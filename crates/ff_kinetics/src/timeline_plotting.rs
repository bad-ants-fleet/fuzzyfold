use std::path::Path;
use rustc_hash::FxHashMap;
use plotters::prelude::*;
use plotters::style::Palette99;

use ff_energy::EnergyModel;
use crate::timeline::Timeline;

pub fn plot_occupancy_over_time<E: EnergyModel>(
    timeline: &Timeline<E>, 
    filename: impl AsRef<Path>,
    title: &str,
    t_lin: f64,
    t_log: f64,
) {
    assert!(t_lin > 0.0 && t_log > t_lin, "Require 0 < t_lin < t_log");

    // Image size; tweak as you like
    // let root = BitMapBackend::new(filename, (1024, 480)).into_drawing_area();
    let root = SVGBackend::new(filename.as_ref(), (1024, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    root.titled(title, ("sans-serif", 28)).unwrap();
    root.draw_text(
        "time",
        &("sans-serif", 22).into_font().into_text_style(&root),
        (496, 450), // roughly centered at bottom
    ).unwrap();


    let eps = 1e-12; // epsilon for plot labels
    // Split into two panels: 50% for linear (left), 50% for log (right)
    let (left, right) = root.split_horizontally(512);

    // ---- Left: linear panel ----
    let mut chart_left = ChartBuilder::on(&left)
        .caption("Linear plot", ("sans-serif", 18))
        .margin(20)
        .margin_top(40)
        .margin_right(0)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..(t_lin+eps), 0.0..1.0).unwrap();
    chart_left
        .configure_mesh()
        .y_desc("occupancy")
        .x_labels(6)
        .x_label_formatter(&|x| if *x < 0.01 {format!("{:.1e}", x)} else {format!("{}", x)})  // scientific notation
        .y_labels(10)
        .light_line_style(RGBColor(220, 220, 220))
        .light_line_style(TRANSPARENT)
        .axis_desc_style(("sans-serif", 22))
        .label_style(("sans-serif", 18))
        .draw()
        .unwrap();

    // draw separator at x = t_lin (right edge of this panel)
    chart_left.draw_series(std::iter::once(PathElement::new(
        vec![(t_lin, 0.0), (t_lin, 1.0)],
        BLACK.mix(0.7),
    ))).unwrap();

    // ---- Right: log panel ----
    let mut chart_right = ChartBuilder::on(&right)
        .caption("Logarithmic plot", ("sans-serif", 18))
        .margin(20)
        .margin_top(40)
        .margin_left(0)
        .margin_right(40)
        .x_label_area_size(40)
        .y_label_area_size(0) // hide y labels on right
        .build_cartesian_2d(((t_lin+eps)..(t_log + eps)).log_scale(), 0.0..1.0)
        .unwrap();

    chart_right
        .configure_mesh()
        .x_labels(6)
        .x_label_formatter(&|x| if *x < 0.01 {format!("{:.1e}", x)} else {format!("{}", x)})  // scientific notation
        .y_labels(10) // hide y ticks on right
        .light_line_style(RGBColor(220, 220, 220))
        .light_line_style(TRANSPARENT)
        .label_style(("sans-serif", 18))
        .draw().unwrap();

    // repeat separator at x = t_lin (left edge of this panel)
    chart_right.draw_series(std::iter::once(PathElement::new(
        vec![(t_lin, 0.0), (t_lin, 1.0)],
        BLACK.mix(0.7),
    ))).unwrap();


    // Group indices by macrostate name, preserving first-seen order.
    let macrostates = timeline.registry.macrostates();
    let mut name_order: Vec<&str> = Vec::new();
    let mut name_to_indices: FxHashMap<&str, Vec<usize>> = FxHashMap::default();
    for (idx, (_len, ms)) in macrostates.iter().enumerate() {
        let name = ms.name();
        name_to_indices
            .entry(name)
            .or_insert_with(|| {
                name_order.push(name);
                Vec::new()
            })
        .push(idx);
        }

    // Build data per macrostate name
    #[allow(clippy::type_complexity)]
    let mut trajectories: Vec<(&str, Vec<(f64, f64, f64)>)> = Vec::new();
    for name in &name_order {
        let indices = &name_to_indices[name];
        let mut series = Vec::new();
        for tp in &timeline.points {
            // Same invariant as the Display impl: at most one length-variant
            // should carry nonzero count at any given timepoint, for now.
            let nonzero_variants = indices.iter()
                .filter(|&&i| tp.ensemble.get(&i).copied().unwrap_or(0) > 0)
                .count();
            debug_assert!(
                nonzero_variants <= 1,
                "macrostate '{}' has nonzero count at {} length-variants \
                simultaneously (indices {:?}) at t={}",
                name, nonzero_variants, indices, tp.time
            );
            let count: usize = indices.iter()
                .map(|&i| tp.ensemble.get(&i).copied().unwrap_or(0))
                .sum();
            let occu = if tp.counter > 0 {
                count as f64 / tp.counter as f64
            } else { 0.0 };
            let se = if tp.counter > 0 {
                (occu * (1.0 - occu) / tp.counter as f64).sqrt()
            } else { 0.0 };
            series.push((tp.time, occu, se));
        }
        if *name == macrostates[0].1.name() || series.iter().any(|(_, occu, _)| *occu >= 0.02) {
            trajectories.push((name, series));
        }
    }

    // Sort by ID to have consistent colors
    //trajectories.sort_by_key(|(id, _)| *id);

    // Find global Y max for normalization
    for (i, (name, series)) in trajectories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9); // pick a distinct color

        let z = 1.0; // or 1.96 for 95%
        let band_color = color.mix(0.2);

        let upper = series.iter().map(|(t, p, se)| (*t, (p + z * se).min(1.0)));
        let lower = series.iter().rev().map(|(t, p, se)| (*t, (p - z * se).max(0.0)));

        let upper = upper.chain(lower);

        chart_left.draw_series(AreaSeries::new(
                upper.clone()
                .filter(|(t, _)| *t <= t_lin),
                0.0,
                band_color,
        )).unwrap();

        chart_left.draw_series(LineSeries::new(
                series.iter().cloned().map(|(t, p, _)| (t, p)).filter(|(t, _)| *t <= t_lin + eps),
                color.stroke_width(2),
        )).unwrap();

        chart_right.draw_series(AreaSeries::new(
                upper
                .filter(|(t, _)| *t >= t_lin),
                0.0,
                band_color,
        )).unwrap();
 
        chart_right.draw_series(LineSeries::new(
            series.iter().cloned().map(|(t, p, _)| (t, p)).filter(|(t, _)| *t >= t_lin - eps),
            color.stroke_width(2),
        )).unwrap()
            .label(name.to_string())   // <-- label for legend
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            });
    }
    
    // after loop:
    chart_right
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
            .label_font(("sans-serif", 16).into_font())   // <-- legend font size
        .draw().unwrap();
    
    root.present().unwrap(); // write the PNG
}
