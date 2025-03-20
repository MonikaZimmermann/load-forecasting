# Mid-Term Hourly Load Forecasting - Zimmermann & Ziel (2025)

This repository contains example R code for the mid-term hourly load forecasting model described in **Zimmermann & Ziel (2025, https://doi.org/10.1016/j.apenergy.2025.125444)**. It demonstrates how to apply the forecasting model using data from **France (FR) and Germany (DE)**. 

The repository provides:
- **Preprocessed data** (load, temperature, holidays, and seasonal patterns) for **24 European countries** covering the period **01/2015–02/2024**, as described in our paper.
- **R scripts** to perform mid-term hourly electricity load forecasting.

If you use this code, please **cite the following paper**:

Zimmermann & Ziel, *"Efficient mid-term forecasting of hourly electricity load using generalized additive models,"* Applied Energy, 2025. DOI: [https://doi.org/10.1016/j.apenergy.2025.125444](https://doi.org/10.1016/j.apenergy.2025.125444)

---

## **Country Code Mapping**
The following table maps the country codes used in the code to their respective country names:

| Code | Country Name       | Code | Country Name       |
|------|--------------------|------|--------------------|
| AT   | Austria            | NL   | Netherlands        |
| BE   | Belgium            | PL   | Poland            |
| BG   | Bulgaria           | PT   | Portugal          |
| CZ   | Czech Republic     | RO   | Romania           |
| DK   | Denmark            | RS   | Serbia            |
| EE   | Estonia            | SK   | Slovakia          |
| FI   | Finland            | SI   | Slovenia          |
| FR   | France             | ES   | Spain             |
| DE   | Germany            | SE   | Sweden            |
| GR   | Greece             | HU   | Hungary           |
| IT   | Italy              | LV   | Latvia            |
| LT   | Lithuania          | ME   | Montenegro        |

---

## **Repository Structure**
- `code/` → Contains R scripts for forecasting functions and example workflows.
- `data/` → Contains preprocessed CSV files with hourly load, temperature, holidays, and seasonal data for 24 European countries.
- `plots/` → Contains generated plots illustrating forecasting results.

---

## **Usage**
To get started, clone the repository


```bash
git clone https://github.com/MonikaZimmermann/load-forecasting.git
cd public-load-forecasting
