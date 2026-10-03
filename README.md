# AEX Conditional Volatility: GARCH vs GJR-GARCH

This project models the time-varying volatility of daily AEX index log-returns from January 2010 to March 2023. It tests whether negative shocks raise volatility more than positive shocks of the same size, which is known as the leverage effect.

**Paper:** [`paper/aex_conditional_volatility.pdf`](paper/aex_conditional_volatility.pdf)

**Authors:** Maxim Milde, Zahid Pashayev (Charles University, 2026)

## Approach

1. **Stationarity:** ADF and KPSS tests on price levels and log-returns.
2. **Mean equation:** ARMA order selection by BIC (`auto.arima` plus manual ARMA(1,1), (2,2) and (5,5)). The result is a constant-mean ARMA(0,0).
3. **ARCH effects:** ACF/PACF of squared residuals, the ARCH-LM test and Ljung-Box on squared residuals.
4. **Volatility models** (`rugarch`): GARCH(1,1) with normal and Student-t errors, and GJR-GARCH(1,1) with Student-t errors.
5. **Model comparison:** Information criteria (AIC, BIC, Shibata, Hannan-Quinn), a likelihood-ratio test for the nested asymmetry term, sign-bias and Nyblom stability tests, and news impact curves.

## Key findings

- Volatility clustering and ARCH effects are strong. Persistence (α + β) is about 0.987.
- Returns are fat-tailed, with excess kurtosis of 6.8. Student-t errors (ν ≈ 7.3) beat normal errors on every information criterion.
- The leverage effect is significant: γ = 0.110, robust t = 2.97. A negative shock moves next-day variance about 6 times as much as a positive shock of the same size (α + γ = 0.132 against α = 0.023).
- GJR-GARCH(1,1)-t is the best-fitting model on all four criteria, and the likelihood-ratio test rejects the symmetric restriction (p < 0.001).

## Reproduce

1. Download daily AEX closing prices (`^AEX`) from Yahoo Finance for 2010-01-04 to 2023-03-06 and save them as `^AEX.csv` in the working directory, with columns `Index` and `AEX.Close`.
2. Run `R/aex_garch_analysis.R`.

Packages used: `tidyverse`, `tseries`, `xts`, `forecast`, `moments`, `aTSA`, `rugarch`.
