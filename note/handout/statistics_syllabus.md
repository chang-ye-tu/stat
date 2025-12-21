# 台大資管所考試統計學準備大綱

**參考書籍：** Jeffrey S. Rosenthal, *Probability and Statistics: The Science of Uncertainty*, 2nd Edition  
**書籍資源：** https://utstat.utoronto.ca/mikevans/jeffrosenthal/

> **說明：** 本大綱根據台大商研、財金、國企、經濟等研究所歷年統計學考古題歸納而成，針對資管所「資管概論」中可能出現的統計學考題進行準備。重點放在線性迴歸模型與假設檢定。

---

## 第零部分：預備微積分知識

> 統計學中的公式推導需要以下微積分基礎，本教材將在各單元中運用這些知識。

### 0.1 基本運算與符號

- **求和符號 (Summation)**
  - $\sum_{i=1}^{n} x_i = x_1 + x_2 + \cdots + x_n$
  - 性質：$\sum (ax_i + b) = a\sum x_i + nb$
  
- **連乘符號 (Product)**
  - $\prod_{i=1}^{n} x_i = x_1 \cdot x_2 \cdots x_n$

- **階乘與組合**
  - $n! = n \cdot (n-1) \cdot (n-2) \cdots 2 \cdot 1$
  - $\binom{n}{k} = \frac{n!}{k!(n-k)!}$

**Rosenthal 參考：** Chapter 1 (Section 1.3)

### 0.2 指數與對數函數

- **指數函數性質**
  - $e^{a+b} = e^a \cdot e^b$
  - $e^{-x} = 1/e^x$
  - $\frac{d}{dx}e^x = e^x$
  - $\frac{d}{dx}e^{ax} = a \cdot e^{ax}$

- **自然對數性質**
  - $\ln(ab) = \ln a + \ln b$
  - $\ln(a/b) = \ln a - \ln b$
  - $\ln(a^n) = n \ln a$
  - $\frac{d}{dx}\ln x = \frac{1}{x}$

- **重要極限**
  - $\lim_{n \to \infty} \left(1 + \frac{x}{n}\right)^n = e^x$
  - $\lim_{x \to 0} \frac{e^x - 1}{x} = 1$

### 0.3 微分基礎

- **基本微分公式**
  - $(x^n)' = nx^{n-1}$
  - $(e^x)' = e^x$
  - $(\ln x)' = 1/x$
  - $(fg)' = f'g + fg'$ (乘法法則)
  - $(f/g)' = (f'g - fg')/g^2$ (除法法則)
  - $(f \circ g)' = f'(g) \cdot g'$ (連鎖法則)

- **偏微分**
  - $\frac{\partial}{\partial x} f(x, y)$：將 $y$ 視為常數對 $x$ 微分
  - 應用於多變量函數的極值問題（如最小平方法推導）

### 0.4 積分基礎

- **不定積分**
  - $\int x^n \, dx = \frac{x^{n+1}}{n+1} + C$ (當 $n \neq -1$)
  - $\int e^{ax} \, dx = \frac{1}{a} e^{ax} + C$
  - $\int \frac{1}{x} \, dx = \ln|x| + C$

- **定積分**
  - $\int_a^b f(x) \, dx = F(b) - F(a)$，其中 $F'(x) = f(x)$

- **換元積分法**
  - 令 $u = g(x)$，則 $\int f(g(x)) g'(x) \, dx = \int f(u) \, du$

- **Gamma 函數**
  - $\Gamma(n) = \int_0^{\infty} x^{n-1} e^{-x} \, dx$
  - $\Gamma(n) = (n-1)!$ 當 $n$ 為正整數
  - $\Gamma(1/2) = \sqrt{\pi}$

**統計應用：** 
- 機率密度函數的積分（求期望值、變異數）
- 常態分配的推導
- 最小平方法的數學證明

---

## 第一部分：機率基礎

### 1.1 機率公理與基本定義

- **樣本空間 (Sample Space)** $\Omega$：所有可能結果的集合
- **事件 (Event)**：樣本空間的子集合

- **機率公理 (Kolmogorov Axioms)**
  1. $0 \leq P(A) \leq 1$
  2. $P(\Omega) = 1$
  3. 若 $A_1, A_2, \ldots$ 互斥，則 $P\left(\bigcup_i A_i\right) = \sum_i P(A_i)$

- **基本性質**
  - $P(A^c) = 1 - P(A)$
  - $P(A \cup B) = P(A) + P(B) - P(A \cap B)$
  - $P(\emptyset) = 0$

**Rosenthal 參考：** Chapter 2 (Sections 2.1–2.3)

### 1.2 條件機率與獨立性

- **條件機率 (Conditional Probability)**
  $$P(A|B) = \frac{P(A \cap B)}{P(B)}, \quad P(B) > 0$$

- **乘法法則**
  $$P(A \cap B) = P(A|B) \cdot P(B) = P(B|A) \cdot P(A)$$

- **獨立性 (Independence)**
  - 事件 $A$ 與 $B$ 獨立 $\Leftrightarrow$ $P(A \cap B) = P(A) \cdot P(B)$
  - 等價於 $P(A|B) = P(A)$

**Rosenthal 參考：** Chapter 2 (Sections 2.4–2.5)

### 1.3 貝氏定理與全機率定理

- **全機率定理 (Law of Total Probability)**
  
  若 $B_1, B_2, \ldots, B_n$ 為樣本空間的分割，則
  $$P(A) = \sum_{i=1}^{n} P(A|B_i) P(B_i)$$

- **貝氏定理 (Bayes' Theorem)**
  $$P(B_j|A) = \frac{P(A|B_j) P(B_j)}{\sum_{i=1}^{n} P(A|B_i) P(B_i)}$$

- **術語**
  - $P(B_j)$：先驗機率 (prior probability)
  - $P(B_j|A)$：後驗機率 (posterior probability)
  - $P(A|B_j)$：似然 (likelihood)

**考古題應用範例：**
- 多選題答對機率問題（知道答案 vs 猜測）
- 購物行為與顧客特徵的條件機率

**Rosenthal 參考：** Chapter 2 (Section 2.6)

### 1.4 聯合機率與邊際機率

- **聯合機率 (Joint Probability)**
  - $P(X = x, Y = y) = P(X = x \cap Y = y)$
  - 聯合機率表所有格子加總等於 1

- **邊際機率 (Marginal Probability)**
  - $P(X = x) = \sum_y P(X = x, Y = y)$
  - 由聯合機率表的行或列加總得到

- **獨立性判定**
  - $X$ 與 $Y$ 獨立 $\Leftrightarrow$ $P(X=x, Y=y) = P(X=x) \cdot P(Y=y)$ 對所有 $x, y$

**Rosenthal 參考：** Chapter 4 (Section 4.1)

---

## 第二部分：隨機變數與常用分配

### 2.1 隨機變數基本概念

- **定義**：隨機變數是從樣本空間到實數的函數
- **離散型 vs 連續型**
  - 離散型：取有限或可數無限個值
  - 連續型：取連續區間內的值

- **機率質量函數 (PMF)** — 離散型
  $$p(x) = P(X = x)$$
  性質：$\sum_x p(x) = 1$

- **機率密度函數 (PDF)** — 連續型
  $$P(a \leq X \leq b) = \int_a^b f(x) \, dx$$
  性質：$\int_{-\infty}^{\infty} f(x) \, dx = 1$，$f(x) \geq 0$

- **累積分配函數 (CDF)**
  $$F(x) = P(X \leq x)$$
  性質：
  - $0 \leq F(x) \leq 1$
  - $F$ 為單調遞增
  - $\lim_{x \to -\infty} F(x) = 0$, $\lim_{x \to \infty} F(x) = 1$

**Rosenthal 參考：** Chapter 3 (Sections 3.1–3.3)

### 2.2 期望值與變異數

- **期望值 (Expectation)**
  - 離散型：$E[X] = \sum_x x \cdot p(x)$
  - 連續型：$E[X] = \int_{-\infty}^{\infty} x \cdot f(x) \, dx$

- **期望值性質**
  - $E[aX + b] = aE[X] + b$
  - $E[X + Y] = E[X] + E[Y]$（無論是否獨立）
  - 若 $X, Y$ 獨立，$E[XY] = E[X] \cdot E[Y]$

- **變異數 (Variance)**
  $$\text{Var}(X) = E[(X - \mu)^2] = E[X^2] - (E[X])^2$$
  其中 $\mu = E[X]$

- **變異數性質**
  - $\text{Var}(aX + b) = a^2 \text{Var}(X)$
  - 若 $X, Y$ 獨立，$\text{Var}(X + Y) = \text{Var}(X) + \text{Var}(Y)$

- **標準差 (Standard Deviation)**
  $$\sigma = \sqrt{\text{Var}(X)}$$

- **共變異數 (Covariance)**
  $$\text{Cov}(X, Y) = E[(X - \mu_X)(Y - \mu_Y)] = E[XY] - E[X]E[Y]$$

- **相關係數 (Correlation)**
  $$\rho_{XY} = \frac{\text{Cov}(X, Y)}{\sigma_X \sigma_Y}, \quad -1 \leq \rho \leq 1$$

**Rosenthal 參考：** Chapter 3 (Sections 3.4–3.5), Chapter 4 (Section 4.2)

### 2.3 常用離散分配

#### 伯努利分配 Bernoulli$(p)$

- **情境**：成功/失敗的單次試驗
- **PMF**：$P(X = 1) = p$, $P(X = 0) = 1-p$
- **期望值與變異數**：$E[X] = p$, $\text{Var}(X) = p(1-p)$

#### 二項分配 Binomial$(n, p)$

- **情境**：$n$ 次獨立伯努利試驗中成功的次數
- **PMF**：
  $$P(X = k) = \binom{n}{k} p^k (1-p)^{n-k}, \quad k = 0, 1, \ldots, n$$
- **期望值與變異數**：$E[X] = np$, $\text{Var}(X) = np(1-p)$
- **常態近似**：當 $np > 5$ 且 $n(1-p) > 5$ 時
  $$X \stackrel{\text{approx}}{\sim} N(np, np(1-p))$$

**考古題應用：**
- 商店達成業績目標的機率
- 連續賭局獲獎機率

#### Poisson 分配 Poisson$(\lambda)$

- **情境**：單位時間/空間內事件發生的次數
- **PMF**：
  $$P(X = k) = \frac{e^{-\lambda} \lambda^k}{k!}, \quad k = 0, 1, 2, \ldots$$
- **期望值與變異數**：$E[X] = \text{Var}(X) = \lambda$
- **與二項分配關係**：當 $n$ 大、$p$ 小、$np = \lambda$ 固定時，Binomial$(n,p) \approx$ Poisson$(\lambda)$

**Rosenthal 參考：** Chapter 5 (Sections 5.1–5.4)

### 2.4 常用連續分配

#### 均勻分配 Uniform$(a, b)$

- **PDF**：
  $$f(x) = \frac{1}{b-a}, \quad a \leq x \leq b$$
- **CDF**：
  $$F(x) = \frac{x-a}{b-a}, \quad a \leq x \leq b$$
- **期望值與變異數**：$E[X] = \frac{a+b}{2}$, $\text{Var}(X) = \frac{(b-a)^2}{12}$

#### 指數分配 Exponential$(\lambda)$

- **情境**：Poisson 過程中事件間的等待時間
- **PDF**：
  $$f(x) = \lambda e^{-\lambda x}, \quad x \geq 0$$
- **CDF**：
  $$F(x) = 1 - e^{-\lambda x}, \quad x \geq 0$$
- **期望值與變異數**：$E[X] = \frac{1}{\lambda}$, $\text{Var}(X) = \frac{1}{\lambda^2}$
- **中位數**：
  $$m = \frac{\ln 2}{\lambda}$$

**推導（考古題重點）：**
$$F(m) = 0.5 \Rightarrow 1 - e^{-\lambda m} = 0.5 \Rightarrow e^{-\lambda m} = 0.5 \Rightarrow m = \frac{\ln 2}{\lambda}$$

- **無記憶性**：$P(X > s + t | X > s) = P(X > t)$

#### 常態分配 Normal$(\mu, \sigma^2)$

- **PDF**：
  $$f(x) = \frac{1}{\sqrt{2\pi}\sigma} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$$
- **期望值與變異數**：$E[X] = \mu$, $\text{Var}(X) = \sigma^2$

- **標準常態分配** $Z \sim N(0, 1)$
  - 標準化：若 $X \sim N(\mu, \sigma^2)$，則 $Z = \frac{X - \mu}{\sigma} \sim N(0, 1)$
  - CDF 記為 $\Phi(z)$

- **常態分配性質**
  - 對稱：$P(Z > z) = P(Z < -z) = 1 - \Phi(z)$
  - 線性組合：若 $X_i \sim N(\mu_i, \sigma_i^2)$ 獨立，則 $\sum a_i X_i \sim N\left(\sum a_i \mu_i, \sum a_i^2 \sigma_i^2\right)$

- **常用分位數（查表）**
  - $z_{0.025} = 1.96$（95% 信賴區間）
  - $z_{0.05} = 1.645$（90% 信賴區間）
  - $z_{0.005} = 2.576$（99% 信賴區間）

**Rosenthal 參考：** Chapter 5 (Sections 5.5–5.8)

---

## 第三部分：抽樣分配

### 3.1 樣本統計量

- **母體 (Population)** vs **樣本 (Sample)**
- **隨機樣本**：$X_1, X_2, \ldots, X_n$ i.i.d.（獨立同分配）

- **常用統計量**
  - 樣本平均數：$\bar{X} = \frac{1}{n}\sum_{i=1}^n X_i$
  - 樣本變異數：$S^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2$
  - 樣本標準差：$S = \sqrt{S^2}$

- **樣本平均數的性質**
  - $E[\bar{X}] = \mu$（不偏估計）
  - $\text{Var}(\bar{X}) = \frac{\sigma^2}{n}$
  - 標準誤 (Standard Error)：$\text{SE}(\bar{X}) = \frac{\sigma}{\sqrt{n}}$

**Rosenthal 參考：** Chapter 6 (Sections 6.1–6.2)

### 3.2 中央極限定理 (Central Limit Theorem)

**定理：** 若 $X_1, X_2, \ldots, X_n$ 為 i.i.d. 隨機變數，$E[X_i] = \mu$，$\text{Var}(X_i) = \sigma^2 < \infty$，則當 $n$ 夠大時：
$$\frac{\bar{X} - \mu}{\sigma/\sqrt{n}} \stackrel{d}{\longrightarrow} N(0, 1)$$

或等價地：
$$\bar{X} \stackrel{\text{approx}}{\sim} N\left(\mu, \frac{\sigma^2}{n}\right)$$

**重要應用：**
- 即使母體非常態，樣本平均數近似常態
- 經驗法則：$n \geq 30$ 通常足夠
- 二項分配的常態近似

**考古題應用範例：**
- IQ 平均值的機率計算
- 樣本平均超過特定值的機率

**Rosenthal 參考：** Chapter 6 (Section 6.3)

### 3.3 抽樣分配

#### 常態母體的抽樣分配

若 $X_1, \ldots, X_n \stackrel{\text{i.i.d.}}{\sim} N(\mu, \sigma^2)$：

1. **樣本平均數分配**
   $$\bar{X} \sim N\left(\mu, \frac{\sigma^2}{n}\right)$$

2. **$\chi^2$ 分配**
   $$\frac{(n-1)S^2}{\sigma^2} \sim \chi^2_{n-1}$$
   - 自由度 $df = n-1$
   - $E[\chi^2_k] = k$, $\text{Var}(\chi^2_k) = 2k$

3. **$t$ 分配**
   $$T = \frac{\bar{X} - \mu}{S/\sqrt{n}} \sim t_{n-1}$$
   - 自由度 $df = n-1$
   - 當 $df$ 大時，$t$ 分配趨近標準常態
   - 常用於母體變異數未知時

4. **$F$ 分配**（兩個獨立 $\chi^2$ 的比）
   $$F = \frac{S_1^2/\sigma_1^2}{S_2^2/\sigma_2^2} \sim F_{n_1-1, n_2-1}$$

**Rosenthal 參考：** Chapter 6 (Sections 6.4–6.6)

---

## 第四部分：點估計與區間估計

### 4.1 點估計

- **估計量 (Estimator)**：用樣本計算的統計量，用來估計母體參數
- **估計值 (Estimate)**：估計量在特定樣本下的值

**良好估計量的性質：**

1. **不偏性 (Unbiasedness)**
   $$E[\hat{\theta}] = \theta$$
   - $\bar{X}$ 是 $\mu$ 的不偏估計量
   - $S^2$ 是 $\sigma^2$ 的不偏估計量

2. **有效性 (Efficiency)**
   - 變異數越小越好
   - 最小變異數不偏估計量 (MVUE)

3. **一致性 (Consistency)**
   $$\hat{\theta}_n \stackrel{p}{\longrightarrow} \theta \text{ as } n \to \infty$$

**Rosenthal 參考：** Chapter 7 (Sections 7.1–7.3)

### 4.2 區間估計

**信賴區間 (Confidence Interval)** 的解釋：
- $(1-\alpha) \times 100\%$ 信賴區間表示：若重複抽樣多次，約有 $(1-\alpha) \times 100\%$ 的區間會包含真正的母體參數

#### 母體平均數的信賴區間

**情況 1：$\sigma$ 已知**
$$\bar{X} \pm z_{\alpha/2} \cdot \frac{\sigma}{\sqrt{n}}$$

**情況 2：$\sigma$ 未知（實務常見）**
$$\bar{X} \pm t_{\alpha/2, n-1} \cdot \frac{S}{\sqrt{n}}$$

#### 母體比例的信賴區間

$$\hat{p} \pm z_{\alpha/2} \cdot \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}$$

其中 $\hat{p} = \frac{X}{n}$（樣本比例）

**誤差界限 (Margin of Error)**
$$E = z_{\alpha/2} \cdot \sqrt{\frac{p(1-p)}{n}}$$

**Rosenthal 參考：** Chapter 7 (Sections 7.4–7.6)

### 4.3 樣本大小決定

#### 估計平均數所需樣本大小

$$n = \left(\frac{z_{\alpha/2} \cdot \sigma}{E}\right)^2$$

其中 $E$ 為所需的誤差界限

#### 估計比例所需樣本大小

$$n = \left(\frac{z_{\alpha/2}}{E}\right)^2 \cdot p(1-p)$$

**保守估計**：使用 $p = 0.5$（使 $p(1-p)$ 最大化）
$$n = \frac{(z_{\alpha/2})^2}{4E^2}$$

**考古題應用範例：**
- 金融公司調查顧客偏好所需樣本數
- 在特定信心水準和誤差範圍下決定樣本大小

**Rosenthal 參考：** Chapter 7 (Section 7.7)

---

## 第五部分：假設檢定

### 5.1 假設檢定基本概念

- **虛無假設 (Null Hypothesis)** $H_0$：欲檢驗的假設，通常包含等號
- **對立假設 (Alternative Hypothesis)** $H_1$：與 $H_0$ 相對的假設

**檢定類型：**
- 雙尾檢定：$H_0: \mu = \mu_0$ vs $H_1: \mu \neq \mu_0$
- 右尾檢定：$H_0: \mu \leq \mu_0$ vs $H_1: \mu > \mu_0$
- 左尾檢定：$H_0: \mu \geq \mu_0$ vs $H_1: \mu < \mu_0$

### 5.2 檢定錯誤與檢定力

|  | $H_0$ 為真 | $H_0$ 為假 |
|---|---|---|
| 不拒絕 $H_0$ | 正確決策 | Type II Error ($\beta$) |
| 拒絕 $H_0$ | Type I Error ($\alpha$) | 正確決策 (Power) |

- **型一錯誤 (Type I Error)**：$H_0$ 為真卻拒絕 $H_0$
  - 機率 = $\alpha$（顯著水準）
  
- **型二錯誤 (Type II Error)**：$H_0$ 為假卻不拒絕 $H_0$
  - 機率 = $\beta$
  
- **檢定力 (Power)**：正確拒絕錯誤 $H_0$ 的機率
  - Power $= 1 - \beta$

**型二錯誤機率計算（考古題重點）：**
1. 在 $H_0$ 下找出接受域
2. 在 $H_1$ 下計算樣本平均落在接受域的機率

**Rosenthal 參考：** Chapter 8 (Sections 8.1–8.3)

### 5.3 p 值

**定義**：在 $H_0$ 為真的假設下，觀察到比當前樣本統計量更極端的機率

**決策規則**：
- 若 p-value $< \alpha$，拒絕 $H_0$
- 若 p-value $\geq \alpha$，不拒絕 $H_0$

**解釋注意事項**：
- p-value 不是 $H_0$ 為真的機率
- p-value 越小，反對 $H_0$ 的證據越強

### 5.4 單樣本檢定

#### $z$ 檢定（$\sigma$ 已知）

**檢定統計量**：
$$z = \frac{\bar{X} - \mu_0}{\sigma/\sqrt{n}}$$

#### $t$ 檢定（$\sigma$ 未知）

**檢定統計量**：
$$t = \frac{\bar{X} - \mu_0}{S/\sqrt{n}}, \quad df = n-1$$

**決策規則（雙尾）**：若 $|t| > t_{\alpha/2, n-1}$，拒絕 $H_0$

**考古題應用範例：**
- 便利商店庫存量是否偏離目標值
- 產品損壞率是否超過標準

**Rosenthal 參考：** Chapter 8 (Section 8.4)

### 5.5 雙樣本檢定

#### 獨立樣本 $t$ 檢定

**假設**：$H_0: \mu_1 - \mu_2 = 0$ vs $H_1: \mu_1 - \mu_2 \neq 0$

**合併變異數估計（假設 $\sigma_1 = \sigma_2$）**：
$$S_p^2 = \frac{(n_1 - 1)S_1^2 + (n_2 - 1)S_2^2}{n_1 + n_2 - 2}$$

**檢定統計量**：
$$t = \frac{(\bar{X}_1 - \bar{X}_2) - 0}{S_p\sqrt{\frac{1}{n_1} + \frac{1}{n_2}}}, \quad df = n_1 + n_2 - 2$$

#### 配對樣本 $t$ 檢定

**情境**：同一受試者的前後測量、配對比較

**步驟**：
1. 計算差異 $d_i = X_{1i} - X_{2i}$
2. 檢定 $H_0: \mu_d = 0$

**檢定統計量**：
$$t = \frac{\bar{d} - 0}{S_d/\sqrt{n}}, \quad df = n - 1$$

**考古題應用範例：**
- 新舊引擎燃油效率比較
- 配對設計 vs 獨立設計的檢定力比較

**Rosenthal 參考：** Chapter 8 (Section 8.5)

### 5.6 比例檢定

**單一比例**：$H_0: p = p_0$

**檢定統計量**：
$$z = \frac{\hat{p} - p_0}{\sqrt{\frac{p_0(1-p_0)}{n}}}$$

**Rosenthal 參考：** Chapter 8 (Section 8.6)

### 5.7 卡方獨立性檢定

**情境**：檢驗兩個類別變數是否獨立

**假設**：
- $H_0$：兩變數獨立
- $H_1$：兩變數不獨立

**期望次數**：
$$E_{ij} = \frac{\text{Row}_i \times \text{Column}_j}{n}$$

**檢定統計量**：
$$\chi^2 = \sum_{i,j} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}$$

**自由度**：$df = (r-1)(c-1)$，其中 $r$ 為列數，$c$ 為行數

**決策**：若 $\chi^2 > \chi^2_{\alpha, df}$，拒絕 $H_0$

**考古題應用範例：**
- 工作級別與年終獎金的關聯性

**Rosenthal 參考：** Chapter 9 (Section 9.3)

---

## 第六部分：簡單線性迴歸

### 6.1 迴歸模型

**模型設定**：
$$Y_i = \beta_0 + \beta_1 X_i + \epsilon_i$$

其中：
- $Y_i$：依變數 (dependent variable)
- $X_i$：自變數 (independent variable)
- $\beta_0$：截距 (intercept)
- $\beta_1$：斜率 (slope)
- $\epsilon_i$：誤差項 (error term)

**迴歸模型假設**：
1. **線性關係**：$E[Y|X] = \beta_0 + \beta_1 X$
2. **誤差期望值為零**：$E[\epsilon_i] = 0$
3. **同質變異數 (Homoscedasticity)**：$\text{Var}(\epsilon_i) = \sigma^2$（常數）
4. **誤差獨立**：$\text{Cov}(\epsilon_i, \epsilon_j) = 0$ for $i \neq j$
5. **常態誤差**（用於假設檢定）：$\epsilon_i \sim N(0, \sigma^2)$

**參數意義**：
- $\beta_0$：當 $X = 0$ 時，$Y$ 的期望值
- $\beta_1$：$X$ 增加 1 單位時，$Y$ 的期望變化量

**Rosenthal 參考：** Chapter 10 (Section 10.1)

### 6.2 最小平方法 (OLS)

**目標**：最小化殘差平方和
$$\text{SSE} = \sum_{i=1}^{n} (Y_i - \hat{Y}_i)^2 = \sum_{i=1}^{n} (Y_i - \hat{\beta}_0 - \hat{\beta}_1 X_i)^2$$

**OLS 估計量的推導**（運用微積分）：

對 $\text{SSE}$ 分別對 $\beta_0$ 和 $\beta_1$ 取偏導數並令其為零：

$$\frac{\partial \text{SSE}}{\partial \beta_0} = -2\sum(Y_i - \beta_0 - \beta_1 X_i) = 0$$
$$\frac{\partial \text{SSE}}{\partial \beta_1} = -2\sum X_i(Y_i - \beta_0 - \beta_1 X_i) = 0$$

**解得**：
$$\hat{\beta}_1 = \frac{\sum_{i=1}^{n}(X_i - \bar{X})(Y_i - \bar{Y})}{\sum_{i=1}^{n}(X_i - \bar{X})^2} = \frac{S_{XY}}{S_{XX}}$$

$$\hat{\beta}_0 = \bar{Y} - \hat{\beta}_1 \bar{X}$$

**與相關係數的關係**：
$$\hat{\beta}_1 = r \cdot \frac{S_Y}{S_X}$$

其中 $r$ 為樣本相關係數，$S_X$ 和 $S_Y$ 為樣本標準差

**考古題重點**：若交換自變數與依變數：
$$\hat{\gamma}_1 = \hat{\beta}_1 \cdot \frac{S_X^2}{S_Y^2}$$

且 $\hat{\beta}_1 \cdot \hat{\gamma}_1 = r^2$

**Rosenthal 參考：** Chapter 10 (Section 10.2)

### 6.3 OLS 估計量的性質

**在模型假設成立時**：

1. **不偏性 (Unbiasedness)**
   $$E[\hat{\beta}_0] = \beta_0, \quad E[\hat{\beta}_1] = \beta_1$$

2. **變異數公式**
   $$\text{Var}(\hat{\beta}_1) = \frac{\sigma^2}{\sum(X_i - \bar{X})^2} = \frac{\sigma^2}{S_{XX}}$$
   $$\text{Var}(\hat{\beta}_0) = \sigma^2 \left(\frac{1}{n} + \frac{\bar{X}^2}{S_{XX}}\right)$$

3. **Gauss-Markov 定理**
   - 在誤差不相關且同質變異數的條件下
   - OLS 是最佳線性不偏估計量 (BLUE: Best Linear Unbiased Estimator)

**重要考古題觀念**：
- 誤差項相關（$\rho \neq 0$）不影響不偏性，但影響效率和變異數估計
- 完全多重共線性會導致 OLS 無法計算

**Rosenthal 參考：** Chapter 10 (Section 10.3)

### 6.4 斜率的假設檢定

**假設**：$H_0: \beta_1 = 0$ vs $H_1: \beta_1 \neq 0$

（$H_0$ 表示 $X$ 與 $Y$ 無線性關係）

**檢定統計量**：
$$t = \frac{\hat{\beta}_1 - 0}{\text{SE}(\hat{\beta}_1)}, \quad df = n - 2$$

其中
$$\text{SE}(\hat{\beta}_1) = \frac{S}{\sqrt{S_{XX}}}, \quad S = \sqrt{\frac{\text{SSE}}{n-2}}$$

**決策**：若 $|t| > t_{\alpha/2, n-2}$，拒絕 $H_0$，表示有顯著線性關係

**Rosenthal 參考：** Chapter 10 (Section 10.4)

### 6.5 斜率的信賴區間

$$\hat{\beta}_1 \pm t_{\alpha/2, n-2} \cdot \text{SE}(\hat{\beta}_1)$$

**變化量的信賴區間**：當 $X$ 變化 $\Delta X$ 時
$$\Delta Y = \hat{\beta}_1 \cdot \Delta X$$

其 $(1-\alpha)$ 信賴區間為：
$$\hat{\beta}_1 \cdot \Delta X \pm t_{\alpha/2, n-2} \cdot |\Delta X| \cdot \text{SE}(\hat{\beta}_1)$$

**考古題應用範例：**
- 高中 GPA 增加對大學 GPA 影響的信賴區間

**Rosenthal 參考：** Chapter 10 (Section 10.4)

### 6.6 判定係數 $R^2$

**定義**：
$$R^2 = \frac{\text{SSR}}{\text{SST}} = 1 - \frac{\text{SSE}}{\text{SST}}$$

其中：
- SST (Total Sum of Squares) $= \sum(Y_i - \bar{Y})^2$
- SSR (Regression Sum of Squares) $= \sum(\hat{Y}_i - \bar{Y})^2$
- SSE (Error Sum of Squares) $= \sum(Y_i - \hat{Y}_i)^2$

**關係式**：SST = SSR + SSE

**性質**：
- $0 \leq R^2 \leq 1$
- $R^2$ 表示 $Y$ 的變異中可被 $X$ 解釋的比例
- 簡單迴歸中：$R^2 = r^2$（相關係數的平方）

**Rosenthal 參考：** Chapter 10 (Section 10.5)

### 6.7 穩健標準誤 (Heteroskedasticity-Robust Standard Errors)

**Eicker-White 估計量**：當誤差變異數不齊一（heteroskedasticity）時使用

$$\hat{D} = \left(\frac{1}{n}\sum \mathbf{x}_i \mathbf{x}_i^\top\right)^{-1} \left(\frac{1}{n}\sum \hat{\epsilon}_i^2 \mathbf{x}_i \mathbf{x}_i^\top\right) \left(\frac{1}{n}\sum \mathbf{x}_i \mathbf{x}_i^\top\right)^{-1}$$

**考古題觀念**：
- 穩健標準誤在同質變異數下仍為一致估計量
- 異質變異數不影響 OLS 的不偏性

---

## 第七部分：迴歸的因果推論議題

### 7.1 相關不等於因果

**關鍵概念**：統計上的顯著關聯不代表因果關係

**需考慮的問題**：
1. **選擇偏誤 (Selection Bias)**：自變數的分配可能與其他變數相關
2. **遺漏變數偏誤 (Omitted Variable Bias)**：重要解釋變數未納入模型
3. **反向因果 (Reverse Causality)**：可能是 Y 影響 X 而非 X 影響 Y
4. **共同原因 (Confounding)**：第三變數同時影響 X 和 Y

### 7.2 識別因果效應的方法

**隨機對照實驗 (RCT)**
- 隨機分配處理組和對照組
- 消除選擇偏誤

**差異中之差異法 (Difference-in-Differences)**
- 需要處理前後的資料
- 需要處理組和對照組
- 模型：$Y = \beta_0 + \beta_1 \cdot \text{Treated} + \beta_2 \cdot \text{Post} + \beta_3 \cdot (\text{Treated} \times \text{Post}) + \epsilon$
- $\beta_3$ 為 DID 估計量

**考古題應用範例（球場案例）**：
- 僅用橫截面資料無法識別因果
- 僅用時間序列資料無法排除時間趨勢
- 需結合時間和橫截面維度的 DID 方法

**Rosenthal 參考：** 補充閱讀（計量經濟學教科書）

---

## 附錄 A：常用統計表摘要

### 標準常態分配 $Z$

| $\alpha$ | $z_\alpha$ |
|----------|-----------|
| 0.10 | 1.282 |
| 0.05 | 1.645 |
| 0.025 | 1.960 |
| 0.01 | 2.326 |
| 0.005 | 2.576 |

### $t$ 分配臨界值

| df | $t_{0.10}$ | $t_{0.05}$ | $t_{0.025}$ | $t_{0.01}$ |
|----|-----------|-----------|------------|-----------|
| 9 | 1.383 | 1.833 | 2.262 | 2.821 |
| 19 | 1.328 | 1.729 | 2.093 | 2.539 |
| 29 | 1.311 | 1.699 | 2.045 | 2.462 |
| 60 | 1.296 | 1.671 | 2.000 | 2.390 |
| $\infty$ | 1.282 | 1.645 | 1.960 | 2.326 |

### $\chi^2$ 分配臨界值

| df | $\chi^2_{0.95}$ | $\chi^2_{0.05}$ | $\chi^2_{0.01}$ |
|----|----------------|----------------|----------------|
| 1 | 0.004 | 3.841 | 6.635 |
| 2 | 0.103 | 5.991 | 9.210 |
| 4 | 0.711 | 9.488 | 13.277 |
| 9 | 3.325 | 16.919 | 21.666 |

---

## 附錄 B：公式速查表

### 機率

- $P(A|B) = P(A \cap B)/P(B)$
- $P(A) = \sum_i P(A|B_i)P(B_i)$
- $P(B_j|A) = P(A|B_j)P(B_j)/P(A)$

### 期望值與變異數

- $E[aX+b] = aE[X] + b$
- $\text{Var}(aX+b) = a^2\text{Var}(X)$
- $\text{Var}(X) = E[X^2] - (E[X])^2$

### 常用分配摘要

| 分配 | 期望值 | 變異數 |
|------|--------|--------|
| Bernoulli$(p)$ | $p$ | $p(1-p)$ |
| Binomial$(n,p)$ | $np$ | $np(1-p)$ |
| Poisson$(\lambda)$ | $\lambda$ | $\lambda$ |
| Uniform$(a,b)$ | $(a+b)/2$ | $(b-a)^2/12$ |
| Exponential$(\lambda)$ | $1/\lambda$ | $1/\lambda^2$ |
| Normal$(\mu,\sigma^2)$ | $\mu$ | $\sigma^2$ |

### 信賴區間

- 平均數（$\sigma$ 未知）：$\bar{X} \pm t_{\alpha/2,n-1} \cdot S/\sqrt{n}$
- 比例：$\hat{p} \pm z_{\alpha/2} \cdot \sqrt{\hat{p}(1-\hat{p})/n}$

### 假設檢定

- $z$-統計量：$z = (\bar{X}-\mu_0)/(\sigma/\sqrt{n})$
- $t$-統計量：$t = (\bar{X}-\mu_0)/(S/\sqrt{n})$
- 配對 $t$-統計量：$t = \bar{d}/(S_d/\sqrt{n})$
- $\chi^2$-統計量：$\chi^2 = \sum (O-E)^2/E$

### 迴歸

- $\hat{\beta}_1 = S_{XY}/S_{XX} = r \cdot S_Y/S_X$
- $\hat{\beta}_0 = \bar{Y} - \hat{\beta}_1\bar{X}$
- $R^2 = 1 - \text{SSE}/\text{SST}$

---

## 附錄 C：Rosenthal 教科書章節對照

| 本大綱章節 | Rosenthal 對應章節 |
|-----------|-------------------|
| 1.1-1.4 機率基礎 | Chapter 2 |
| 2.1-2.2 隨機變數與期望值 | Chapter 3 |
| 2.3-2.4 常用分配 | Chapter 5 |
| 3.1-3.3 抽樣分配 | Chapter 6 |
| 4.1-4.3 點估計與區間估計 | Chapter 7 |
| 5.1-5.7 假設檢定 | Chapter 8-9 |
| 6.1-6.7 線性迴歸 | Chapter 10 |

**補充資源：**
- 書籍全文：https://utstat.utoronto.ca/mikevans/jeffrosenthal/book.pdf
- 習題解答：https://utstat.utoronto.ca/mikevans/jeffrosenthal/solutions.pdf

---

## 附錄 D：練習題指引

每個單元建議練習題目：

### 機率基礎
- Rosenthal Ch.2: #1, 3, 5, 7, 11, 15, 19

### 隨機變數
- Rosenthal Ch.3: #1, 3, 5, 9, 13
- Rosenthal Ch.5: #1, 5, 7, 11, 15

### 抽樣分配
- Rosenthal Ch.6: #1, 3, 5, 9

### 估計與檢定
- Rosenthal Ch.7: #1, 3, 5, 7
- Rosenthal Ch.8: #1, 3, 5, 7, 9

### 迴歸分析
- Rosenthal Ch.10: #1, 3, 5, 7

---

*本大綱根據台大研究所歷年統計學考古題編製，持續更新中。*
