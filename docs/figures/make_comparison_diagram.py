# Side-by-side comparison diagram: scMultiomeGRN vs MeVD-GRN (SVG; render with headless Chrome).
import textwrap
BLUE, GREEN, GREY, INK = "#0072B2", "#009E73", "#6B7280", "#1F2937"
W, H = 1600, 860
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Poppins, Inter, Arial, sans-serif">',
     '<defs><filter id="sh" x="-10%" y="-10%" width="120%" height="130%"><feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="#000" flood-opacity="0.10"/></filter></defs>',
     f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>',
     f'<text x="60" y="64" font-size="36" font-weight="700" fill="{INK}">scMultiomeGRN vs MeVD-GRN</text>',
     f'<text x="60" y="98" font-size="18" fill="{GREY}">Same ingredients (paired RNA + ATAC, graph encoders, link-prediction loss), different design choices</text>']
cols = ["TASK", "INPUT GRAPH", "ENCODER", "MODALITY FUSION", "DECODER", "TRAINING LABELS"]
A = ["Undirected TF-TF link completion", "The known (training) edges", "GraFRank attention + edge histograms", "Cross-modal attention, per node", "Concat-MLP (paper says inner product)", "Motif hits from the same ATAC (circular)"]
B = ["Directed TF → any gene", "Label-free: co-expression + TF-candidates", "GraphSAGE towers + Geneformer prior", "None until the decoder", "Bilinear TF→gene score gated by ATAC openness", "ChIP, then knockout; their overlap held out"]
x0, cw, gap = 225, 205, 12
for j, c in enumerate(cols):
    o.append(f'<text x="{x0 + j*(cw+gap) + cw/2}" y="150" text-anchor="middle" font-size="15" font-weight="600" letter-spacing="1.4" fill="{GREY}">{c}</text>')
def card(x, y, w, h, color, txt, fs=17):
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="#FFFFFF" stroke="{color}" stroke-width="3" filter="url(#sh)"/>')
    lines = textwrap.wrap(txt, 20); y0 = y + h/2 - (len(lines)-1)*11 + 6
    for k, l in enumerate(lines): o.append(f'<text x="{x+w/2}" y="{y0+k*23}" text-anchor="middle" font-size="{fs}" fill="{INK}">{l}</text>')
for r, (name, color, row) in enumerate([("scMultiomeGRN", GREEN, A), ("MeVD-GRN", BLUE, B)]):
    y = 175 + r*190
    o.append(f'<rect x="60" y="{y}" width="150" height="150" rx="16" fill="{color}"/>')
    o.append(f'<text x="135" y="{y+80}" text-anchor="middle" font-size="17" font-weight="700" fill="#FFFFFF">{name}</text>')
    for j, t in enumerate(row): card(x0 + j*(cw+gap), y, cw, 150, color, t)
# training strip
o.append(f'<text x="60" y="590" font-size="15" font-weight="600" letter-spacing="1.4" fill="{GREY}">TRAINING</text>')
for r, (name, color, txt) in enumerate([("scMultiomeGRN", GREEN, "10 runs, up to 2,000 epochs each; keep edges that win a 6-of-10 vote; one flat label set"),
                                       ("MeVD-GRN", BLUE, "45 epochs, no voting; evidence curriculum (ChIP → knockout); zero-shot test on the ChIP ∩ knockout tier")]):
    y = 610 + r*90
    o.append(f'<rect x="60" y="{y}" width="150" height="70" rx="16" fill="{color}"/><text x="135" y="{y+42}" text-anchor="middle" font-size="15" font-weight="700" fill="#FFFFFF">{name}</text>')
    o.append(f'<rect x="{x0}" y="{y}" width="{6*cw+5*gap}" height="70" rx="16" fill="#FFFFFF" stroke="{color}" stroke-width="3" filter="url(#sh)"/>')
    o.append(f'<text x="{x0+22}" y="{y+42}" font-size="18" fill="{INK}">{txt}</text>')
o.append(f'<text x="60" y="815" font-size="15" fill="{GREY}">Shared: GraFRank-style separate RNA and ATAC pathways; class-weighted BCE loss; MAESTRO-style regulatory-potential ATAC features.</text>')
o.append('</svg>')
open("docs/figures/comparison_scmultiomegrn_vs_mevd.svg", "w").write("\n".join(o)); print("ok")
