use plotters::prelude::*;
use crate::timeline::Macrostate;
use plotters::style::Palette99;
use crate::timeline::Timeline;

pub fn plot_occupancy_over_time(
    timeline: &Timeline, 
    filename: &str,
    t_lin: f64,
    t_log: f64,
) {
    assert!(t_lin > 0.0 && t_log > t_lin, "Require 0 < t_lin < t_log");

    let title = format!("ff-simulate ({} simulations)", timeline.points[0].counter);

    // Image size; tweak as you like
    //let root = BitMapBackend::new(filename, (1024, 480)).into_drawing_area();
    let root = SVGBackend::new(filename, (1024, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    root.titled(&title, ("sans-serif", 28)).unwrap();
    root.draw_text(
        "time",
        &("sans-serif", 18).into_font().into_text_style(&root),
        (512, 450),   // roughly centered at bottom
    ).unwrap();


    // Split into two panels: 18% for linear (left), 82% for log (right)
    let (left, right) = root.split_horizontally(200);

    // ---- Left: linear panel ----
    let mut chart_left = ChartBuilder::on(&left)
        //.caption(format!("Linear plot"), ("sans-serif", 14))
        .margin(20)
        .margin_top(40)
        .margin_right(0)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..t_lin, 0.0..1.0).unwrap();
    chart_left
        .configure_mesh()
        //.x_desc("time (s)")
        .y_desc("occupancy")
        .x_labels(2)
        .x_label_formatter(&|x| {
            if (*x - t_lin).abs() < 1e-9 { "".to_string() } else { format!("{}", x) }
        })
    .y_labels(10)
        .light_line_style(&RGBColor(220, 220, 220))
        .axis_desc_style(("sans-serif", 18))
        .label_style(("sans-serif", 18))
        .draw()
        .unwrap();

    // draw separator at x = t_lin (right edge of this panel)
    chart_left.draw_series(std::iter::once(PathElement::new(
        vec![(t_lin, 0.0), (t_lin, 1.0)],
        &BLACK.mix(0.7),
    ))).unwrap();

    // ---- Right: log panel ----
    let mut chart_right = ChartBuilder::on(&right)
        //.caption(format!("Logarithmic plot"), ("sans-serif", 14))
        .margin(20)
        .margin_top(40)
        .margin_left(0)
        .margin_right(40)
        .x_label_area_size(40)
        .y_label_area_size(0) // hide y labels on right
        .build_cartesian_2d((t_lin..t_log).log_scale(), 0.0..1.0).unwrap();

    chart_right
        .configure_mesh()
        //.x_desc("time (s)")
        //.x_labels(7)
        .x_label_formatter(&|x| if *x < 0.01 {format!("{:.1e}", x)} else {format!("{}", x)})  // scientific notation
        .y_labels(10) // hide y ticks on right
        .light_line_style(&RGBColor(220, 220, 220))
        .label_style(("sans-serif", 18))
        .draw().unwrap();

    // repeat separator at x = t_lin (left edge of this panel)
    chart_right.draw_series(std::iter::once(PathElement::new(
        vec![(t_lin, 0.0), (t_lin, 1.0)],
        &BLACK.mix(0.7),
    ))).unwrap();


    // Build data per structure
    let mut trajectories: Vec<(usize, Vec<(f64, f64)>)> = Vec::new();

    for (id, _) in timeline.registry.iter() {
        let mut series = Vec::new();
        for tp in &timeline.points {
            let count = tp.ensemble.get(&id).copied().unwrap_or(0);
            let occu = if tp.counter > 0 {
                count as f64 / tp.counter as f64
            } else {
                0.0
            };
            series.push((tp.time, occu));
        }
        if series.iter().any(|(_, occu)| *occu >= 0.1) { // threshold filter
            trajectories.push((id, series));
        }
    }

    // Sort by ID to have consistent colors
    trajectories.sort_by_key(|(id, _)| *id);

    // Find global Y max for normalization
    for (i, (id, series)) in trajectories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9); // pick a distinct color

        let (name, energy) = match timeline.registry.get(*id) {
            Macrostate::Explicit(reg) => (reg.name(), reg.energy()),
            Macrostate::Constraint(_) => ("Constraint", None),
        };

        chart_left.draw_series(LineSeries::new(
                series.iter().cloned().filter(|(t, _)| *t <= t_lin),
                color.stroke_width(2),
        )).unwrap();

        chart_right.draw_series(LineSeries::new(
            series.iter().cloned().filter(|(t, _)| *t >= t_lin),
            color.stroke_width(2),
        )).unwrap()
            .label(format!("{:20} {:>6.2}", name.trim(), energy.unwrap_or(0.0)))   // <-- label for legend
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            });
    }
    
// after loop:
chart_right
    .configure_series_labels()
    .border_style(&BLACK)
    .background_style(&WHITE.mix(0.8))
    .position(SeriesLabelPosition::UpperRight)
        .label_font(("sans-serif", 16).into_font())   // <-- legend font size
    .draw().unwrap();

    root.present().unwrap(); // write the PNG
}


