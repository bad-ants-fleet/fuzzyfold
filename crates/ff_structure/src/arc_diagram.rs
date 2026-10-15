use crate::PairTable;

pub fn arc_diagram_svg(sequence: &str, pt: &PairTable) -> String {
    arc_diagram_svg_opt(Some(sequence), pt)
}

/// Arc diagram with an optional sequence. Without one, bases are drawn as
/// neutral dots and no nucleotide labels are emitted.
pub fn arc_diagram_svg_opt(sequence: Option<&str>, pt: &PairTable) -> String {
    let n = pt.len();
    if let Some(seq) = sequence {
        assert_eq!(seq.len(), n, "sequence length must match pair table length");
    }

    let margin = 40.0_f64;
    let base_spacing = 20.0_f64;
    let base_r = 7.0_f64;
    let font_size = 11.0_f64;

    let x_pos = |i: usize| margin + i as f64 * base_spacing;

    // --- collect arcs (only i < j) ---
    let mut arcs: Vec<(usize, usize)> = Vec::new();
    for (i, entry) in pt.iter().enumerate() {
        if let Some(j) = entry {
            let j = *j as usize;
            if i < j {
                arcs.push((i, j));
            }
        }
    }

    // Natural height for each arc: span/2 = semicircle, never looks flat
    let natural_heights: Vec<f64> = arcs
        .iter()
        .map(|&(i, j)| {
            let xi = x_pos(i);
            let xj = x_pos(j);
            (xj - xi) / 2.0
        })
        .collect();

    // base_y must accommodate the tallest arc + margin
    let max_natural_h = natural_heights.iter().cloned().fold(0.0_f64, f64::max);
    let base_y = max_natural_h + margin;

    let height = base_y + margin + base_r + font_size + 4.0;
    let width = margin * 2.0 + (n as f64 - 1.0) * base_spacing;

    // Rank arcs by span for the small disambiguation offset
    let disambiguation_step = 4.0_f64; // px per rank step
    let mut ranked: Vec<usize> = (0..arcs.len()).collect();
    ranked.sort_by_key(|&idx| {
        let (i, j) = arcs[idx];
        j - i
    });
    let mut arc_rank = vec![0usize; arcs.len()];
    for (rank, &idx) in ranked.iter().enumerate() {
        arc_rank[idx] = rank;
    }

    // Final height = natural (semicircle) + small rank offset
    let arc_heights: Vec<f64> = (0..arcs.len())
        .map(|idx| natural_heights[idx] + arc_rank[idx] as f64 * disambiguation_step)
        .collect();

    // --- SVG ---
    let mut svg = String::new();

    // viewBox lets consumers scale the diagram; the explicit rect is the actual
    // background -- CSS `background` on the root is not honoured by rasterizers.
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="monospace"><rect width="{width}" height="{height}" fill="#ffffff"/>"##,
        width = width, height = height
    ));

    // Backbone
    svg.push_str(&format!(
        r##"<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" stroke="#585b70" stroke-width="1.5"/>"##,
        x1 = margin,
        y = base_y,
        x2 = margin + (n as f64 - 1.0) * base_spacing,
    ));

    // Arcs
    for (idx, (i, j)) in arcs.iter().enumerate() {
        let xi = x_pos(*i);
        let xj = x_pos(*j);
        let cx = (xi + xj) / 2.0;
        let arc_h = arc_heights[idx];

        svg.push_str(&format!(
            r##"<path d="M {xi} {y} Q {cx} {apex} {xj} {y}" fill="none" stroke="#89b4fa" stroke-width="1.6" stroke-opacity="0.85"/>"##,
            xi = xi, xj = xj, y = base_y,
            cx = cx, apex = base_y - arc_h,
        ));
    }

    // Base circles + labels
    let colors = [
        ("A", "#a6e3a1"),
        ("U", "#f38ba8"),
        ("G", "#fab387"),
        ("C", "#89dceb"),
    ];
    let base_color = |ch: char| -> &'static str {
        colors
            .iter()
            .find(|(b, _)| b.starts_with(ch))
            .map(|(_, c)| *c)
            .unwrap_or("#cdd6f4")
    };

    for i in 0..n {
        let ch = sequence.and_then(|s| s.chars().nth(i));
        let x = x_pos(i);
        let color = match ch {
            Some(c) => base_color(c.to_ascii_uppercase()),
            None => "#cdd6f4",
        };
        let stroke = if pt[i].is_some() {
            "#ffffff"
        } else {
            "#585b70"
        };

        svg.push_str(&format!(
            r##"<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>"##,
            x = x, y = base_y, r = if ch.is_some() { base_r } else { base_r * 0.6 },
            fill = color, stroke = stroke,
        ));
        if let Some(c) = ch {
            svg.push_str(&format!(
                r##"<text x="{x}" y="{ty}" text-anchor="middle" font-size="{fs}" fill="#313244">{ch}</text>"##,
                x = x, ty = base_y + base_r + font_size + 2.0, fs = font_size, ch = c,
            ));
        }
        if i % 5 == 0 {
            svg.push_str(&format!(
                r##"<text x="{x}" y="{ty}" text-anchor="middle" font-size="9" fill="#6c7086">{i}</text>"##,
                x = x, ty = base_y - base_r - 4.0, i = i,
            ));
        }
    }

    svg.push_str("</svg>");
    svg
}

/// Generates the arc diagram and writes it to a file.
pub fn save_arc_diagram(
    sequence: &str,
    pt: &PairTable,
    output_filename: impl AsRef<std::path::Path>,
) -> Result<(), std::io::Error> {
    let svg = arc_diagram_svg(sequence, pt);
    std::fs::write(output_filename, svg)
}

/// As [`save_arc_diagram`], with an optional sequence.
pub fn save_arc_diagram_opt(
    sequence: Option<&str>,
    pt: &PairTable,
    output_filename: impl AsRef<std::path::Path>,
) -> Result<(), std::io::Error> {
    std::fs::write(output_filename, arc_diagram_svg_opt(sequence, pt))
}
