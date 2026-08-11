APP_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=DM+Sans:wght@300;400;500&display=swap');

html, body, [class*="css"] { font-family: 'DM Sans', sans-serif; }
code, .stCode, pre { font-family: 'DM Mono', monospace !important; }

/* Top bar */
.block-container { padding-top: 2rem; padding-bottom: 2rem; }

/* Metric cards */
[data-testid="metric-container"] {
    background: #f8f7f4;
    border: 1px solid #e8e5de;
    border-radius: 10px;
    padding: 1rem;
}
[data-testid="stMetricLabel"] { font-size: 0.72rem; color: #888; letter-spacing: 0.06em; text-transform: uppercase; }
[data-testid="stMetricValue"] { font-size: 1.5rem; font-weight: 500; }

/* Table */
[data-testid="stDataFrame"] { border: 1px solid #e8e5de; border-radius: 10px; overflow: hidden; }

/* Upload box */
[data-testid="stFileUploadDropzone"] {
    border: 1.5px dashed #c8c5be !important;
    border-radius: 10px !important;
    background: #fafaf8 !important;
}

/* Buttons */
.stDownloadButton > button {
    background: #1a1a1a !important;
    color: #fff !important;
    border: none !important;
    border-radius: 8px !important;
    font-size: 0.85rem !important;
    padding: 0.5rem 1.2rem !important;
    font-family: 'DM Sans', sans-serif !important;
}
.stDownloadButton > button:hover { opacity: 0.85 !important; }

/* Status badges */
.badge-stable   { background:#e6f4ec; color:#1a7a3e; padding:2px 10px; border-radius:20px; font-size:0.78rem; }
.badge-unstable { background:#fdecea; color:#b91c1c; padding:2px 10px; border-radius:20px; font-size:0.78rem; }

h1 { font-weight: 500 !important; letter-spacing: -0.02em; }
h3 { font-weight: 500 !important; color: #444; }
</style>
"""

EXPECTED_FORMAT_EXAMPLE = (
    "sequence,mutation_region,num_random_copies\n"
    "MPYEKHVEQTVVEKTE,3-5,4\n"
    "RQGGGAPAGGNIGGG,,1\n"
    "MKTQRDGHSLGRWSLV,6,8"
)
