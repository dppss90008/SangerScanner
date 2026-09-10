# SangerScanner

一個桌面 GUI 工具，用來讀取 Sanger 定序的 **AB1 (.ab1)** 原始檔，
自動偵測**雙峰（heterozygous / diploid）訊號位置**、做 **IUPAC 混合鹼基**轉換，
並提供互動式的訊號軌跡、雙峰位置清單與定序品質圖。

適用情境：單一模板定序時想快速找出可能的雜合位點、確認雙峰、輸出 IUPAC consensus 序列
（例如念珠菌 MLST、真菌 / 細菌鑑定、SNP 檢查等）。

---

## 功能特色

- **讀取 AB1 原始檔**，取出四色（A/C/G/T）訊號軌跡與 basecaller 資訊
- **雙峰偵測**：以「次峰振幅 ÷ 主峰振幅」比值 ≥ 門檻（預設 `0.25`）判定為雙峰位置
- **IUPAC 轉換**：雙峰位置自動轉成混合鹼基代碼（R/Y/S/W/K/M/B/D/H/V/N）
- **三個互動分頁**
  - **序列訊號** — 四色軌跡疊圖，可調整每段 bp、拖曳捲動、右側顯示 IUPAC 序列（可切換反向互補 RC）
  - **雙峰訊號位置** — 次峰比例散點圖 + 雙峰位置清單，點圖或點列會同步標示，右側顯示該位置局部放大圖
  - **定序品質** — 每個位置的品質分數（Phred-like，取自 ABIF `PCON2`）+ 雙峰位置清單
- **一鍵分析並自動匯出**：按「開始分析」後，會在 `.ab1` 檔旁自動建立一個以檔名命名的資料夾，
  把所有結果（FASTA、CSV、圖檔）存進去
- **反向互補（Reverse Complement）**：正確處理 IUPAC 混合鹼基，並一併匯出
- 文字（檔案路徑、狀態列、序列面板）皆可選取 / 右鍵複製

---

## 安裝與執行

### 方式一：直接用 Windows exe（推薦給一般使用者）

不需要安裝 Python。

1. 到本專案的 **[Actions](../../actions)** 分頁，點最新一次成功（綠勾）的 *Build Windows exe*
2. 下方 **Artifacts** 下載 `SangerScanner-windows`（zip）
3. 解壓縮得到 **`SangerScanner.exe`**，雙擊即可執行

### 方式二：從原始碼執行（Mac / Linux / Windows 皆可）

需要 Python 3.9–3.12。

```bash
pip install -r requirements.txt
python ab1_gui.py
```

（`requirements.txt` 內含 `pyinstaller`，只執行原始碼的話可以不裝它。）

### 方式三：自己在 Windows 上打包 exe

把整個資料夾複製到 Windows 電腦，雙擊 **`build_exe.bat`**。
它會建立一個乾淨的虛擬環境、安裝相依套件、用 PyInstaller 打包，
完成後 exe 位於 `dist\SangerScanner.exe`。

---

## 使用方式

1. 按 **瀏覽…** 選擇一個 `.ab1` 檔
2. 需要的話調整 **雙峰訊號門檻**（預設 0.25，數字越小越敏感、假陽性越多）
3. 按 **開始分析**
4. 在三個分頁檢視結果；結果檔案已同時自動匯出到 `.ab1` 檔旁的同名資料夾

> 提醒：Sanger 讀長的**開頭與結尾**訊噪比通常較差，容易出現假的雙峰。
> 請搭配「定序品質」分頁與局部放大圖，肉眼確認每個雙峰位置是否可信。

---

## 分析方法

| 項目 | 說明 |
|---|---|
| 訊號軌跡 | ABIF `DATA9–DATA12`（已做 baseline / mobility 校正），依 `FWO_1` 對應到 A/C/G/T |
| 主要序列 | basecaller 判定的 `PBAS2` |
| 主峰 / 次峰振幅 | `P1AM1` / `P2AM1`；次要鹼基 `P2BA1` |
| 雙峰判定 | `peak_ratio = P2AM1 / P1AM1`，`peak_ratio ≥ 門檻` → 判為雙峰 |
| IUPAC 轉換 | 雙峰位置 → 主要+次要鹼基合併成 IUPAC 代碼；其餘位置維持原鹼基 |
| 品質分數 | ABIF `PCON2`，Phred-like（`Q = -10·log10(P_error)`），由定序儀 basecaller 產生 |

---

## 輸出檔案

分析後會在 `.ab1` 檔所在目錄建立 `<檔名>/` 資料夾，內含：

| 檔案 | 內容 |
|---|---|
| `<檔名>_IUPAC.fasta` | IUPAC consensus 序列（雙峰位置為混合鹼基） |
| `<檔名>_IUPAC_RC.fasta` | 上者的反向互補序列 |
| `<檔名>_het_calls.csv` | 每個位置的主/次鹼基、主/次峰振幅、比值、品質、是否雙峰、IUPAC |
| `<檔名>_full_trace_stacked.png` | 整條序列的四色軌跡疊圖（分段堆疊成一張長圖） |
| `<檔名>_peak_ratio_overview.png` | 次峰比例總覽散點圖 |
| `<檔名>_quality_overview.png` | 每位置品質分數圖 |
| `<檔名>_het_zoom.png` | 每個雙峰位置的局部放大圖集 |

---

## 專案檔案

| 檔案 | 用途 |
|---|---|
| `ab1_gui.py` | 主程式（單一檔） |
| `requirements.txt` | Python 相依套件 |
| `build_exe.bat` | 在 Windows 上打包 exe |
| `.github/workflows/build-windows.yml` | GitHub Actions 雲端打包 Windows exe |
| `C26 9-56_CPAR2_808110-F.ab1` | 測試用範例 AB1 檔 |

---

## 相依套件

`biopython`、`numpy`、`pandas`、`matplotlib`（`tkinter` 為 Python 內建）。

---

## 作者

Author: CH Hsieh — SangerScanner v1.0.0
