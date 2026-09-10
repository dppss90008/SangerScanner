#!/usr/bin/env python3
"""
SangerScanner - Sanger AB1 Diploid / IUPAC 判讀 GUI

用法:
    python3 ab1_gui.py

需要套件: biopython, numpy, pandas, matplotlib (tkinter 為 Python 內建)
本機已驗證可用的直譯器: /Users/ch.hsieh/anaconda3/bin/python3
"""
import os
import re
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("缺少 biopython，請先執行: pip install biopython")


BASE_COLORS = {"A": "#2E7D32", "C": "#1565C0", "G": "#424242", "T": "#C62828"}

IUPAC_CODE = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}

matplotlib.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
    # 讓圖上的中文提示文字（例如放大圖面板的預設提示）能正確顯示，不要缺字
    "font.sans-serif": ["Microsoft JhengHei", "Microsoft YaHei", "PingFang TC",
                        "Heiti TC", "Arial Unicode MS", "DejaVu Sans"],
    "axes.unicode_minus": False,
})


def iupac_call(b1, b2):
    if b1 == b2:
        return b1
    return IUPAC_CODE[frozenset([b1, b2])]


IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def reverse_complement(seq):
    """回傳 IUPAC 序列的反向互補 (Reverse Complement)，正確處理混合鹼基代碼。"""
    return "".join(IUPAC_COMPLEMENT.get(b, b) for b in reversed(seq))


_SEQ_LINE_PREFIX_RE = re.compile(r"(?m)^\s*\d+ ")


def strip_seq_line_prefix(text):
    """複製 IUPAC 序列面板內容時，把每行開頭的位置編號去掉，只留下序列本身。"""
    return _SEQ_LINE_PREFIX_RE.sub("", text)


class Ab1Data:
    """讀取單一 ab1 檔案需要的原始欄位。"""

    def __init__(self, path):
        record = SeqIO.read(path, "abi")
        ann = record.annotations["abif_raw"]
        base_order = ann["FWO_1"].decode()

        self.path = path
        self.sample_id = record.id
        self.traces = {
            base_order[0]: np.array(ann["DATA9"], dtype=float),
            base_order[1]: np.array(ann["DATA10"], dtype=float),
            base_order[2]: np.array(ann["DATA11"], dtype=float),
            base_order[3]: np.array(ann["DATA12"], dtype=float),
        }
        self.peak_loc = np.array(ann["PLOC2"])
        self.primary_base = ann["PBAS2"].decode()
        self.secondary_base = ann["P2BA1"].decode()
        self.primary_amp = np.array(ann["P1AM1"], dtype=float)
        self.secondary_amp = np.array(ann["P2AM1"], dtype=float)
        self.quality = np.array(list(ann["PCON2"]))


def analyze(data, het_threshold):
    ratio = np.divide(data.secondary_amp, data.primary_amp,
                       out=np.zeros_like(data.primary_amp), where=data.primary_amp != 0)
    is_het = ratio >= het_threshold

    df = pd.DataFrame({
        "base_index": np.arange(len(data.primary_base)),
        "scan_position": data.peak_loc,
        "primary_base": list(data.primary_base),
        "secondary_base": list(data.secondary_base),
        "primary_amp": data.primary_amp,
        "secondary_amp": data.secondary_amp,
        "peak_ratio": ratio,
        "quality": data.quality,
        "is_het": is_het,
    })
    df["iupac_call"] = [
        iupac_call(r.primary_base, r.secondary_base) if r.is_het else r.primary_base
        for r in df.itertuples()
    ]
    consensus = "".join(df["iupac_call"])
    return df, consensus


def plot_region(data, df, base_start, base_end, ax, pad=15, show_legend=True,
                 label_fontsize=12, axis_fontsize=None, title_fontsize=None, tick_fontsize=None):
    """畫出 base_start~base_end (含) 範圍的四色 trace 疊圖到既有的 ax 上。

    label_fontsize/axis_fontsize/title_fontsize/tick_fontsize 可用來放大字體
    （例如匯出 het-zoom 圖集時），不指定則沿用預設大小。
    """
    ax.clear()
    b0, b1 = max(0, base_start), min(len(df) - 1, base_end)
    scan_start = max(0, int(data.peak_loc[b0]) - pad)
    scan_end = min(len(data.traces["A"]), int(data.peak_loc[b1]) + pad)
    x = np.arange(scan_start, scan_end)

    for base, color in BASE_COLORS.items():
        ax.plot(x, data.traces[base][scan_start:scan_end], color=color, lw=1.5, label=base)

    ymax = max(data.traces[base][scan_start:scan_end].max() for base in "ACGT") * 1.15
    ymax = max(ymax, 10)
    for i in range(b0, b1 + 1):
        row = df.iloc[i]
        label_color = "#C62828" if row.is_het else "#333333"
        weight = "bold" if row.is_het else "normal"
        ax.axvline(row.scan_position, color="#BBBBBB", lw=0.6, linestyle="--", zorder=0)
        ax.text(row.scan_position, ymax, row.iupac_call, ha="center", va="bottom",
                fontsize=label_fontsize, color=label_color, fontweight=weight)

    ax.set_ylim(0, ymax * 1.25)
    ax.set_xlim(scan_start, scan_end)
    ax.set_xlabel("Scan position", fontsize=axis_fontsize)
    ax.set_ylabel("Signal intensity", fontsize=axis_fontsize)
    if tick_fontsize is not None:
        ax.tick_params(labelsize=tick_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if b0 == b1:
        row = df.iloc[b0]
        ax.set_title(f"#{b0}: {row.primary_base}/{row.secondary_base} → {row.iupac_call}  "
                     f"(ratio={row.peak_ratio:.2f}, Q={row.quality})", fontsize=title_fontsize)
    else:
        ax.set_title(f"Bases {b0}-{b1}", fontsize=title_fontsize)

    if show_legend:
        ax.legend(loc="upper right", ncol=4, frameon=False, fontsize=axis_fontsize)


def plot_ratio_overview(df, het_threshold, ax):
    ax.clear()
    het_mask = df["is_het"].values
    ax.scatter(df.loc[~het_mask, "base_index"], df.loc[~het_mask, "peak_ratio"],
               c="#90A4AE", s=18, label="homozygous")
    ax.scatter(df.loc[het_mask, "base_index"], df.loc[het_mask, "peak_ratio"],
               c="#FB8C00", s=26, label="heterozygous (≥ threshold)")
    ax.axhline(het_threshold, color="#424242", lw=1.2, linestyle="--",
               label=f"threshold = {het_threshold}")
    ax.set_xlabel("Base position")
    ax.set_ylabel("Secondary / primary peak ratio")
    ax.set_title("Secondary-peak signal strength across the read")
    ymax = min(df["peak_ratio"].max() * 1.1, 2.0) if len(df) else 1.0
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)


def plot_quality_overview(df, ax):
    ax.clear()
    het_mask = df["is_het"].values
    ax.plot(df["base_index"], df["quality"], color="#CFD8DC", lw=1, zorder=1)
    ax.scatter(df.loc[~het_mask, "base_index"], df.loc[~het_mask, "quality"],
               c="#90A4AE", s=14, label="homozygous", zorder=2)
    ax.scatter(df.loc[het_mask, "base_index"], df.loc[het_mask, "quality"],
               c="#FB8C00", s=22, label="heterozygous", zorder=3)
    ax.axhline(20, color="#424242", lw=1, linestyle="--", label="Q = 20")
    ax.set_xlabel("Base position")
    ax.set_ylabel("Quality (Phred-like, from PCON2)")
    ax.set_title("Per-base quality across the read")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)


class Ab1App(tk.Tk):
    SEQ_PREFIX_LEN = 5   # 每行開頭的位置編號寬度 (含一個空格)

    def __init__(self):
        super().__init__()
        self.title("SangerScanner")
        # 依實際螢幕大小決定視窗尺寸並置中，避免視窗比螢幕大而被切掉右邊/下面的內容
        self.update_idletasks()
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()
        win_w = min(1180, screen_w - 80)
        win_h = min(820, screen_h - 100)
        x = (screen_w - win_w) // 2
        y = max(0, (screen_h - win_h) // 3)
        self.geometry(f"{win_w}x{win_h}+{x}+{y}")
        self.minsize(900, 650)

        self.data = None
        self.df = None
        self.consensus = ""
        self.consensus_rc = ""
        self.seq_view_mode = "fwd"  # "fwd" 或 "rc"，控制序列面板顯示方向
        self.window_size = tk.IntVar(value=60)
        self.het_threshold = tk.DoubleVar(value=0.25)
        self.seq_row_width = tk.IntVar(value=20)  # IUPAC 序列每行顯示幾個 base，可於面板調整
        self.current_start = 0

        self._build_menu()
        self._build_title_bar()
        self._build_top_bar()
        self._build_info_bar()
        self._build_notebook()
        self._build_statusbar()
        self._bind_keys()

    # ---------- UI 建置 ----------
    def _build_menu(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="開啟 AB1 檔案...", command=self.browse_file)
        filemenu.add_separator()
        filemenu.add_command(label="結束", command=self.quit)
        menubar.add_cascade(label="檔案", menu=filemenu)
        self.config(menu=menubar)

    def _build_title_bar(self):
        ttk.Label(self, text="Author: CH Hsieh | SangerScanner-ver1.0.0",
                  foreground="#616161", anchor="center",
                  font=("Menlo", 11)).pack(side=tk.TOP, fill=tk.X, pady=(6, 0))

    def _build_top_bar(self):
        bar = ttk.Frame(self, padding=8)
        bar.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(bar, text="AB1 檔案:").pack(side=tk.LEFT)
        self.path_var = tk.StringVar()
        ttk.Entry(bar, textvariable=self.path_var, width=55).pack(side=tk.LEFT, padx=4)
        ttk.Button(bar, text="瀏覽...", command=self.browse_file).pack(side=tk.LEFT, padx=2)

        ttk.Label(bar, text="   雙峰訊號門檻:").pack(side=tk.LEFT)
        ttk.Entry(bar, textvariable=self.het_threshold, width=6).pack(side=tk.LEFT)

        self.analyze_btn = ttk.Button(bar, text="開始分析", command=self.run_analysis)
        self.analyze_btn.pack(side=tk.RIGHT, padx=2)
        self.progress = ttk.Progressbar(bar, mode="indeterminate", length=100)
        # 分析中才 pack 出來，平常不佔位置

    def _build_info_bar(self):
        info_bar = ttk.Frame(self, padding=(8, 0, 8, 6))
        info_bar.pack(side=tk.TOP, fill=tk.X)
        self.info_var = tk.StringVar(value="尚未分析")
        # 用唯讀 Entry 而非 Label，讓文字可以選取 / 複製 (Cmd/Ctrl-C 或右鍵複製)
        info_entry = ttk.Entry(info_bar, textvariable=self.info_var, state="readonly")
        info_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self._add_copy_context_menu(info_entry)

    def _build_het_treeview(self, parent):
        """建立一個顯示雙峰位置清單的 Treeview（供多個分頁共用同一套欄位設定）。"""
        columns = ("base_index", "primary", "secondary", "primary_amp", "secondary_amp",
                   "iupac", "ratio", "quality")
        headers = {"base_index": "位置", "primary": "主要", "secondary": "次要",
                   "primary_amp": "主要強度", "secondary_amp": "次要強度",
                   "iupac": "IUPAC", "ratio": "比例", "quality": "品質"}
        tree = ttk.Treeview(parent, columns=columns, show="headings")
        for c in columns:
            tree.heading(c, text=headers[c])
            tree.column(c, width=80, anchor="center")
        tscroll = ttk.Scrollbar(parent, orient="vertical", command=tree.yview)
        txscroll = ttk.Scrollbar(parent, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=tscroll.set, xscrollcommand=txscroll.set)
        tscroll.pack(side=tk.RIGHT, fill=tk.Y)
        txscroll.pack(side=tk.BOTTOM, fill=tk.X)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree.bind("<Double-1>", self.on_tree_double_click)
        tree.bind("<<TreeviewSelect>>", self._on_tree_select)
        return tree

    def _build_notebook(self):
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        # -- 區段疊圖瀏覽 --
        self.trace_frame = ttk.Frame(self.nb)
        self.nb.add(self.trace_frame, text="序列訊號")
        paned = ttk.PanedWindow(self.trace_frame, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        plot_container = ttk.Frame(paned)
        paned.add(plot_container, weight=3)

        self.fig_trace = Figure(figsize=(10, 4.6), dpi=100)
        self.ax_trace = self.fig_trace.add_subplot(111)
        toolbar_frame = ttk.Frame(plot_container)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        self.canvas_trace = FigureCanvasTkAgg(self.fig_trace, master=plot_container)
        NavigationToolbar2Tk(self.canvas_trace, toolbar_frame)
        self.canvas_trace.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        nav_bar = ttk.Frame(plot_container)
        nav_bar.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 4))
        ttk.Button(nav_bar, text="◀ 上一段", command=self.prev_region).pack(side=tk.LEFT, padx=4)
        self.region_label = ttk.Label(nav_bar, text="尚未載入")
        self.region_label.pack(side=tk.LEFT, padx=10)
        ttk.Button(nav_bar, text="下一段 ▶", command=self.next_region).pack(side=tk.LEFT, padx=4)

        ttk.Label(nav_bar, text="   每段 bp:").pack(side=tk.LEFT)
        window_spin = ttk.Spinbox(nav_bar, from_=5, to=500, increment=5, width=5,
                                   textvariable=self.window_size,
                                   command=self._on_window_size_change)
        window_spin.pack(side=tk.LEFT)
        window_spin.bind("<Return>", lambda e: self._on_window_size_change())
        window_spin.bind("<FocusOut>", lambda e: self._on_window_size_change())

        scroll_bar = ttk.Frame(plot_container)
        scroll_bar.pack(side=tk.BOTTOM, fill=tk.X, padx=4, pady=(0, 4))
        ttk.Label(scroll_bar, text="拖曳捲動:").pack(side=tk.LEFT)
        self.scale_trace = ttk.Scale(scroll_bar, from_=0, to=1, orient=tk.HORIZONTAL,
                                      command=self._on_scale_change)
        self.scale_trace.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)

        # -- 右側：IUPAC 序列面板 (可與左側疊圖一起拖曳調整寬度) --
        seq_container = ttk.Frame(paned)
        paned.add(seq_container, weight=1)
        seq_container.rowconfigure(1, weight=1)
        seq_container.columnconfigure(0, weight=1)

        seq_header = ttk.Frame(seq_container)
        seq_header.grid(row=0, column=0, sticky="ew", padx=4, pady=(4, 2))
        self.seq_header_var = tk.StringVar(value="IUPAC 序列")
        ttk.Label(seq_header, textvariable=self.seq_header_var).pack(side=tk.LEFT)

        seq_body = ttk.Frame(seq_container)
        seq_body.grid(row=1, column=0, sticky="nsew", padx=4)
        self.seq_text = tk.Text(seq_body, wrap="char", font=("Menlo", 13),
                                 width=self.seq_row_width.get() + self.SEQ_PREFIX_LEN + 2)
        seq_yscroll = ttk.Scrollbar(seq_body, orient="vertical", command=self.seq_text.yview)
        self.seq_text.configure(yscrollcommand=seq_yscroll.set)
        self.seq_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        seq_yscroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.seq_text.tag_configure("het", background="#FFA726", foreground="#000000",
                                     font=("Menlo", 13, "bold"))
        self.seq_text.tag_configure("current", background="#FFF176")
        self.seq_text.tag_configure("posnum", foreground="#9E9E9E")
        self.seq_text.tag_raise("het")  # 雙峰的橘色標示優先於目前範圍的黃色標示
        # 唯讀但仍可選取 / 複製：攔截按鍵輸入，不影響滑鼠選取與 Cmd/Ctrl-C
        self.seq_text.bind("<Key>", lambda e: "break")
        self._add_copy_context_menu(self.seq_text, cleaner=strip_seq_line_prefix)
        # 攔截 Cmd/Ctrl-C，複製時把行首的位置編號去掉，只留序列本身
        self.seq_text.bind("<<Copy>>", self._on_seq_text_copy)

        seq_footer = ttk.Frame(seq_container)
        seq_footer.grid(row=2, column=0, sticky="ew", padx=4, pady=4)
        ttk.Label(seq_footer, text="每行 bases:").pack(side=tk.LEFT)
        seq_width_spin = ttk.Spinbox(seq_footer, from_=5, to=200, increment=5, width=4,
                                      textvariable=self.seq_row_width,
                                      command=self._on_seq_width_change)
        seq_width_spin.pack(side=tk.LEFT, padx=(4, 10))
        seq_width_spin.bind("<Return>", lambda e: self._on_seq_width_change())
        seq_width_spin.bind("<FocusOut>", lambda e: self._on_seq_width_change())
        self.rc_toggle_btn = ttk.Button(seq_footer, text="顯示反向互補 (RC)",
                                         command=self.toggle_rc_view)
        self.rc_toggle_btn.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(seq_footer, text="複製全部序列到剪貼簿", command=self.copy_sequence).pack(
            side=tk.LEFT, fill=tk.X, expand=True)

        # -- 雙峰訊號位置 (合併一頁；圖表與清單互相點選連動) --
        self.table_frame = ttk.Frame(self.nb)
        self.nb.add(self.table_frame, text="雙峰訊號位置")

        paned2 = ttk.PanedWindow(self.table_frame, orient=tk.VERTICAL)
        paned2.pack(fill=tk.BOTH, expand=True)

        chart_container = ttk.Frame(paned2)
        paned2.add(chart_container, weight=3)
        self.fig_overview = Figure(figsize=(10, 4.2), dpi=100)
        self.ax_overview = self.fig_overview.add_subplot(111)
        toolbar_frame2 = ttk.Frame(chart_container)
        toolbar_frame2.pack(side=tk.TOP, fill=tk.X)
        self.canvas_overview = FigureCanvasTkAgg(self.fig_overview, master=chart_container)
        NavigationToolbar2Tk(self.canvas_overview, toolbar_frame2)
        self.canvas_overview.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.overview_annot = None
        self.overview_selected_marker = None
        self.canvas_overview.mpl_connect("motion_notify_event", self._on_overview_hover)
        self.canvas_overview.mpl_connect("button_press_event", self._on_overview_click)

        table_container = ttk.Frame(paned2)
        paned2.add(table_container, weight=2)
        ttk.Label(table_container, padding=(4, 4),
                  text="雙峰位置清單 — 點圖上的點或點列，右邊會顯示該位置的局部放大圖").pack(
            side=tk.TOP, anchor="w")

        table_split = ttk.PanedWindow(table_container, orient=tk.HORIZONTAL)
        table_split.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        table_body = ttk.Frame(table_split)
        table_split.add(table_body, weight=2)
        self.tree = self._build_het_treeview(table_body)

        zoom_container = ttk.Frame(table_split)
        table_split.add(zoom_container, weight=1)
        self.fig_zoom = Figure(figsize=(5, 4), dpi=100)
        self.ax_zoom = self.fig_zoom.add_subplot(111)
        self.ax_zoom.axis("off")
        self.ax_zoom.text(0.5, 0.5, "點選左側或上方圖表\n查看該位置放大圖",
                           ha="center", va="center", color="#9E9E9E", fontsize=11,
                           transform=self.ax_zoom.transAxes)
        self.canvas_zoom = FigureCanvasTkAgg(self.fig_zoom, master=zoom_container)
        self.canvas_zoom.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # -- 定序品質 --
        self.quality_frame = ttk.Frame(self.nb)
        self.nb.add(self.quality_frame, text="定序品質")

        paned3 = ttk.PanedWindow(self.quality_frame, orient=tk.VERTICAL)
        paned3.pack(fill=tk.BOTH, expand=True)

        quality_chart_container = ttk.Frame(paned3)
        paned3.add(quality_chart_container, weight=3)
        self.fig_quality = Figure(figsize=(10, 4.2), dpi=100)
        self.ax_quality = self.fig_quality.add_subplot(111)
        toolbar_q = ttk.Frame(quality_chart_container)
        toolbar_q.pack(side=tk.TOP, fill=tk.X)
        self.canvas_quality = FigureCanvasTkAgg(self.fig_quality, master=quality_chart_container)
        NavigationToolbar2Tk(self.canvas_quality, toolbar_q)
        self.canvas_quality.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.quality_annot = None
        self.quality_selected_marker = None
        self.canvas_quality.mpl_connect("motion_notify_event", self._on_quality_hover)

        quality_table_container = ttk.Frame(paned3)
        paned3.add(quality_table_container, weight=2)
        ttk.Label(quality_table_container, padding=(4, 4),
                  text="雙峰位置清單").pack(side=tk.TOP, anchor="w")
        self.tree_quality = self._build_het_treeview(quality_table_container)

    def _build_statusbar(self):
        self.status_var = tk.StringVar(value="請選擇 AB1 檔案")
        # 用唯讀 Entry 而非 Label，讓文字（例如匯出路徑）可以選取 / 複製 (Cmd/Ctrl-C 或右鍵複製)
        status_entry = ttk.Entry(self, textvariable=self.status_var, state="readonly")
        status_entry.pack(side=tk.BOTTOM, fill=tk.X)
        self._add_copy_context_menu(status_entry)

    def _bind_keys(self):
        self.bind("<Left>", lambda e: self.prev_region())
        self.bind("<Right>", lambda e: self.next_region())

    def _add_copy_context_menu(self, widget, cleaner=None):
        """右鍵選單 (Windows/Linux 用 Button-3，macOS 觸控板/滑鼠右鍵多半也是 Button-2)。
        cleaner: 可選，複製前用來過濾/轉換文字的函式（例如去掉序列面板的行首位置編號）。
        """
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(label="複製", command=lambda: self._copy_widget_content(widget, cleaner))

        def show_menu(event):
            try:
                menu.tk_popup(event.x_root, event.y_root)
            finally:
                menu.grab_release()

        widget.bind("<Button-2>", show_menu)
        widget.bind("<Button-3>", show_menu)

    def _copy_widget_content(self, widget, cleaner=None):
        try:
            text = widget.selection_get()
        except tk.TclError:
            if isinstance(widget, tk.Text):
                text = widget.get("1.0", "end-1c")
            else:
                text = widget.get()
        if cleaner is not None:
            text = cleaner(text)
        self.clipboard_clear()
        self.clipboard_append(text)

    # ---------- 動作 ----------
    def browse_file(self):
        path = filedialog.askopenfilename(
            title="選擇 AB1 檔案", filetypes=[("AB1 files", "*.ab1"), ("All files", "*.*")])
        if path:
            self.path_var.set(path)
            self.run_analysis()

    def run_analysis(self):
        path = self.path_var.get().strip()
        if not path or not os.path.exists(path):
            messagebox.showerror("錯誤", "請選擇有效的 AB1 檔案")
            return
        try:
            threshold = float(self.het_threshold.get())
            window = int(self.window_size.get())
            if window <= 0:
                raise ValueError
        except (tk.TclError, ValueError):
            messagebox.showerror("錯誤", "雙峰訊號門檻需為有效數字")
            return

        self.analyze_btn.config(state="disabled", text="分析中...")
        self.progress.pack(side=tk.RIGHT, padx=(0, 8))
        self.progress.start(10)
        self.status_var.set("讀取並分析中...")
        self.update_idletasks()
        try:
            try:
                self.data = Ab1Data(path)
                self.df, self.consensus = analyze(self.data, threshold)
                self.consensus_rc = reverse_complement(self.consensus)
                self.seq_view_mode = "fwd"
            except Exception as exc:
                messagebox.showerror("讀取失敗", str(exc))
                self.status_var.set("讀取失敗")
                return

            self.current_start = 0
            max_start = max(0, len(self.df) - 1)
            self.scale_trace.configure(to=max_start)
            self.scale_trace.set(0)
            self._refresh_info_bar()
            self._refresh_table()
            self._refresh_overview()
            self._refresh_quality()
            self._populate_sequence_panel()
            self._render_current_region()

            first_iid = self.tree.get_children()
            if first_iid:
                self.tree.selection_set(first_iid[0])
                self.tree.see(first_iid[0])

            het_n = int(self.df["is_het"].sum())

            self.status_var.set("匯出結果中...")
            self.update_idletasks()
            out_dir = self.export_all(silent=True)
            if out_dir is not None:
                self.status_var.set(
                    f"完成: {self.data.sample_id}，共 {len(self.df)} bp，偵測到 {het_n} 個雙峰位置，"
                    f"結果已自動匯出至: {out_dir}")
            else:
                self.status_var.set(
                    f"完成: {self.data.sample_id}，共 {len(self.df)} bp，偵測到 {het_n} 個雙峰位置"
                    f"（自動匯出失敗）")
        finally:
            self.progress.stop()
            self.progress.pack_forget()
            self.analyze_btn.config(state="normal", text="開始分析")

    def _refresh_info_bar(self):
        het_n = int(self.df["is_het"].sum())
        self.info_var.set(
            f"檔案: {self.data.path}    "
            f"雙峰訊號門檻: {self.het_threshold.get()}    "
            f"偵測到雙峰位置: {het_n} 個 (共 {len(self.df)} bp)"
        )

    def _refresh_table(self):
        het_df = self.df[self.df["is_het"]]
        for tree in (self.tree, self.tree_quality):
            tree.delete(*tree.get_children())
            for r in het_df.itertuples():
                tree.insert("", tk.END, iid=str(r.base_index), values=(
                    r.base_index, r.primary_base, r.secondary_base,
                    f"{r.primary_amp:.0f}", f"{r.secondary_amp:.0f}", r.iupac_call,
                    f"{r.peak_ratio:.2f}", r.quality,
                ))

    def _refresh_overview(self):
        plot_ratio_overview(self.df, float(self.het_threshold.get()), self.ax_overview)
        self.overview_annot = self.ax_overview.annotate(
            "", xy=(0, 0), xycoords="data",
            xytext=(15, 15), textcoords="figure pixels",
            bbox=dict(boxstyle="round", fc="#FFF3E0", ec="#FB8C00", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="#FB8C00", lw=1.2),
            fontsize=11, visible=False, zorder=10)
        self.overview_selected_marker = None
        self.fig_overview.tight_layout()
        self.canvas_overview.draw()

        self.ax_zoom.clear()
        self.ax_zoom.axis("off")
        self.ax_zoom.text(0.5, 0.5, "點選左側或上方圖表\n查看該位置放大圖",
                           ha="center", va="center", color="#9E9E9E", fontsize=11,
                           transform=self.ax_zoom.transAxes)
        self.canvas_zoom.draw()

    def _on_overview_hover(self, event):
        if self.overview_annot is None or self.df is None:
            return
        if event.inaxes != self.ax_overview or event.xdata is None:
            if self.overview_annot.get_visible():
                self.overview_annot.set_visible(False)
                self.canvas_overview.draw_idle()
            return
        idx = int(round(event.xdata))
        idx = max(0, min(idx, len(self.df) - 1))
        row = self.df.iloc[idx]
        self.overview_annot.xy = (row.base_index, row.peak_ratio)  # 箭頭指向實際資料點
        self.overview_annot.xyann = (event.x + 15, event.y + 15)   # 文字框貼著滑鼠游標
        self.overview_annot.set_text(
            f"pos {row.base_index}: {row.primary_base}/{row.secondary_base}\n"
            f"ratio={row.peak_ratio:.2f}  Q={row.quality}\n"
            f"IUPAC={row.iupac_call}")
        self.overview_annot.set_visible(True)
        self.canvas_overview.draw_idle()

    def _refresh_quality(self):
        plot_quality_overview(self.df, self.ax_quality)
        self.quality_annot = self.ax_quality.annotate(
            "", xy=(0, 0), xycoords="data",
            xytext=(15, 15), textcoords="figure pixels",
            bbox=dict(boxstyle="round", fc="#FFF3E0", ec="#FB8C00", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="#FB8C00", lw=1.2),
            fontsize=11, visible=False, zorder=10)
        self.quality_selected_marker = None
        self.fig_quality.tight_layout()
        self.canvas_quality.draw()

    def _on_quality_hover(self, event):
        if self.quality_annot is None or self.df is None:
            return
        if event.inaxes != self.ax_quality or event.xdata is None:
            if self.quality_annot.get_visible():
                self.quality_annot.set_visible(False)
                self.canvas_quality.draw_idle()
            return
        idx = int(round(event.xdata))
        idx = max(0, min(idx, len(self.df) - 1))
        row = self.df.iloc[idx]
        self.quality_annot.xy = (row.base_index, row.quality)
        self.quality_annot.xyann = (event.x + 15, event.y + 15)
        self.quality_annot.set_text(
            f"pos {row.base_index}: {row.primary_base}\n"
            f"quality={row.quality}\n"
            f"IUPAC={row.iupac_call}")
        self.quality_annot.set_visible(True)
        self.canvas_quality.draw_idle()

    def _on_overview_click(self, event):
        if self.df is None or event.inaxes != self.ax_overview or event.xdata is None:
            return
        idx = int(round(event.xdata))
        idx = max(0, min(idx, len(self.df) - 1))
        row = self.df.iloc[idx]
        if row.is_het:
            iid = str(idx)
            if self.tree.exists(iid):
                self.tree.selection_set(iid)
                self.tree.see(iid)
        else:
            self.tree.selection_remove(*self.tree.selection())
            self.status_var.set(f"位置 {idx} 非雙峰位置 ({row.primary_base})")

    def _on_tree_select(self, event):
        widget = event.widget
        sel = widget.selection()
        if not sel or self.df is None:
            return
        idx = int(sel[0])
        row = self.df.iloc[idx]

        # <<TreeviewSelect>> 是非同步（排入事件佇列）觸發，用旗標防重入沒有用；
        # 必須先確認「目標樹的選取真的不同」才呼叫 selection_set，
        # 否則兩個 Treeview 會永遠互相觸發對方，造成無窮迴圈 (CPU 100% 當掉)。
        other = self.tree_quality if widget is self.tree else self.tree
        iid = str(idx)
        if other.exists(iid) and other.selection() != (iid,):
            other.selection_set(iid)
            other.see(iid)

        if self.overview_selected_marker is not None:
            try:
                self.overview_selected_marker.remove()
            except ValueError:
                pass
        self.overview_selected_marker = self.ax_overview.scatter(
            [row.base_index], [row.peak_ratio], s=140, facecolors="none",
            edgecolors="#E65100", linewidths=2, zorder=5)
        self.canvas_overview.draw_idle()

        if self.quality_selected_marker is not None:
            try:
                self.quality_selected_marker.remove()
            except ValueError:
                pass
        self.quality_selected_marker = self.ax_quality.scatter(
            [row.base_index], [row.quality], s=140, facecolors="none",
            edgecolors="#E65100", linewidths=2, zorder=5)
        self.canvas_quality.draw_idle()

        plot_region(self.data, self.df, idx, idx, self.ax_zoom, pad=8, show_legend=True)
        self.fig_zoom.tight_layout()
        self.canvas_zoom.draw_idle()

    def _populate_sequence_panel(self):
        self.seq_text.delete("1.0", tk.END)
        width = max(1, int(self.seq_row_width.get()))
        if self.seq_view_mode == "rc":
            seq = self.consensus_rc
            het_mask = self.df["is_het"].values[::-1]
        else:
            seq = self.consensus
            het_mask = self.df["is_het"].values
        for i in range(0, len(seq), width):
            chunk = seq[i:i + width]
            mask_chunk = het_mask[i:i + width]
            self.seq_text.insert(tk.END, f"{i:>4} ", "posnum")
            for ch, is_h in zip(chunk, mask_chunk):
                self.seq_text.insert(tk.END, ch, "het" if is_h else "")
            self.seq_text.insert(tk.END, "\n")

    def _highlight_range(self, b0, b1):
        if self.df is None:
            return
        self.seq_text.tag_remove("current", "1.0", tk.END)
        if self.seq_view_mode == "rc":
            n = len(self.df)
            b0, b1 = n - 1 - b1, n - 1 - b0  # 正向範圍映射到 RC 座標
        width = max(1, int(self.seq_row_width.get()))
        prefix = self.SEQ_PREFIX_LEN
        first_line, last_line = b0 // width, b1 // width
        for line in range(first_line, last_line + 1):
            line_base0 = line * width
            seg0 = max(b0, line_base0)
            seg1 = min(b1, line_base0 + width - 1)
            row = line + 1
            col0 = prefix + (seg0 - line_base0)
            col1 = prefix + (seg1 - line_base0) + 1
            self.seq_text.tag_add("current", f"{row}.{col0}", f"{row}.{col1}")
        self.seq_text.see(f"{first_line + 1}.0")

    def _on_window_size_change(self):
        if self.df is None:
            return
        try:
            window = int(self.window_size.get())
            if window <= 0:
                raise ValueError
        except (tk.TclError, ValueError):
            return
        self._render_current_region()

    def _on_seq_width_change(self):
        if self.df is None:
            return
        try:
            width = int(self.seq_row_width.get())
            if width <= 0:
                raise ValueError
        except (tk.TclError, ValueError):
            return
        self._populate_sequence_panel()
        self._render_current_region()

    def _on_seq_text_copy(self, _event):
        self._copy_widget_content(self.seq_text, cleaner=strip_seq_line_prefix)
        return "break"

    def copy_sequence(self):
        if not self.consensus:
            return
        seq = self.consensus_rc if self.seq_view_mode == "rc" else self.consensus
        self.clipboard_clear()
        self.clipboard_append(seq)
        label = "反向互補 (RC)" if self.seq_view_mode == "rc" else "IUPAC"
        self.status_var.set(f"已複製{label}序列到剪貼簿")

    def toggle_rc_view(self):
        if self.df is None:
            return
        self.seq_view_mode = "rc" if self.seq_view_mode == "fwd" else "fwd"
        if self.seq_view_mode == "rc":
            self.seq_header_var.set("反向互補序列 RC")
            self.rc_toggle_btn.config(text="顯示正向序列 (Forward)")
        else:
            self.seq_header_var.set("IUPAC 序列")
            self.rc_toggle_btn.config(text="顯示反向互補 (RC)")
        self._populate_sequence_panel()
        self._render_current_region()

    def _render_current_region(self):
        if self.df is None:
            return
        window = int(self.window_size.get())
        b0 = self.current_start
        b1 = min(b0 + window - 1, len(self.df) - 1)
        plot_region(self.data, self.df, b0, b1, self.ax_trace)
        self.fig_trace.tight_layout()
        self.canvas_trace.draw()
        self.region_label.config(text=f"Base {b0}-{b1}  (共 {len(self.df)} bp)")
        self._highlight_range(b0, b1)

    def _on_scale_change(self, value_str):
        if self.df is None:
            return
        start = int(round(float(value_str)))
        if start != self.current_start:
            self.current_start = start
            self._render_current_region()

    def _goto(self, start):
        if self.df is None:
            return
        start = max(0, min(start, len(self.df) - 1))
        self.scale_trace.set(start)  # 觸發 _on_scale_change 重繪
        if start == self.current_start:
            self._render_current_region()  # 邊界值未變時仍強制刷新一次

    def next_region(self):
        self._goto(self.current_start + int(self.window_size.get()))

    def prev_region(self):
        self._goto(self.current_start - int(self.window_size.get()))

    def on_tree_double_click(self, event):
        sel = event.widget.selection()
        if not sel or self.df is None:
            return
        base_index = int(sel[0])
        window = int(self.window_size.get())
        self._goto(max(0, base_index - window // 2))
        self.nb.select(self.trace_frame)

    def export_all(self, silent=False):
        if self.df is None:
            if not silent:
                messagebox.showwarning("尚未分析", "請先載入並分析 AB1 檔案")
            return None
        base_name = os.path.splitext(os.path.basename(self.data.path))[0]
        out_dir = os.path.join(os.path.dirname(self.data.path), base_name)
        os.makedirs(out_dir, exist_ok=True)

        out_prefix = os.path.join(out_dir, base_name)
        window = int(self.window_size.get())
        threshold = float(self.het_threshold.get())

        if not silent:
            self.status_var.set("匯出中...")
            self.update_idletasks()
        try:
            fasta_path = f"{out_prefix}_IUPAC.fasta"
            with open(fasta_path, "w") as fh:
                fh.write(f">{self.data.sample_id} diploid_consensus threshold={threshold}\n")
                for i in range(0, len(self.consensus), 60):
                    fh.write(self.consensus[i:i + 60] + "\n")

            rc_fasta_path = f"{out_prefix}_IUPAC_RC.fasta"
            with open(rc_fasta_path, "w") as fh:
                fh.write(f">{self.data.sample_id} diploid_consensus_reverse_complement "
                         f"threshold={threshold}\n")
                for i in range(0, len(self.consensus_rc), 60):
                    fh.write(self.consensus_rc[i:i + 60] + "\n")

            csv_path = f"{out_prefix}_het_calls.csv"
            self.df.to_csv(csv_path, index=False)

            window_starts = list(range(0, len(self.df), window))

            stacked_fig = Figure(figsize=(14, 4.2 * len(window_starts)), dpi=110)

            for row, b0 in enumerate(window_starts):
                b1 = min(b0 + window - 1, len(self.df) - 1)
                ax_stack = stacked_fig.add_subplot(len(window_starts), 1, row + 1)
                plot_region(self.data, self.df, b0, b1, ax_stack, show_legend=(row == 0))

            stacked_fig.tight_layout()
            stacked_fig.savefig(f"{out_prefix}_full_trace_stacked.png", bbox_inches="tight")

            fig_ov = Figure(figsize=(15, 5), dpi=150)
            ax_ov = fig_ov.add_subplot(111)
            plot_ratio_overview(self.df, threshold, ax_ov)
            fig_ov.tight_layout()
            fig_ov.savefig(f"{out_prefix}_peak_ratio_overview.png", bbox_inches="tight")

            fig_q = Figure(figsize=(15, 5), dpi=150)
            ax_q = fig_q.add_subplot(111)
            plot_quality_overview(self.df, ax_q)
            fig_q.tight_layout()
            fig_q.savefig(f"{out_prefix}_quality_overview.png", bbox_inches="tight")

            het_positions = self.df.loc[self.df["is_het"], "base_index"].tolist()
            if het_positions:
                ncols = 3
                nrows = int(np.ceil(len(het_positions) / ncols))
                fig_gallery = Figure(figsize=(6.5 * ncols, 5 * nrows), dpi=110)
                for idx, b in enumerate(het_positions):
                    ax_g = fig_gallery.add_subplot(nrows, ncols, idx + 1)
                    plot_region(self.data, self.df, b, b, ax_g, pad=8, show_legend=False,
                                label_fontsize=20, axis_fontsize=16, title_fontsize=17,
                                tick_fontsize=14)
                fig_gallery.suptitle("Zoom-in view of every heterozygous (IUPAC) call", fontsize=22)
                fig_gallery.tight_layout()
                fig_gallery.savefig(f"{out_prefix}_het_zoom.png", bbox_inches="tight")
        except Exception as exc:
            messagebox.showerror("匯出失敗", str(exc))
            self.status_var.set("匯出失敗")
            return None

        if not silent:
            self.status_var.set(f"已匯出全部結果到: {out_dir}")
            messagebox.showinfo("完成", f"已匯出全部結果到:\n{out_dir}")
        return out_dir


def main():
    app = Ab1App()
    app.mainloop()


if __name__ == "__main__":
    main()
