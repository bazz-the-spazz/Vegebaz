# Calculation of indicator-values and more

This is a collection of R-scripts that:
- ... extract the indicator_value data from
    - [Vegedaz 2021](https://www.wsl.ch/de/services-produkte/vegedaz/) (*read_data_from_Vegedaz.r*)
      - based on Flora Indicativa (Flora indicativa und NISM Mai 2019 Nuller korrigiert)
      - based on Landolt (Zeigerwerte Landolt und NISM Mai 2019)
    - [Vegedaz 2023 Beta Dez](https://www.wsl.ch/de/services-produkte/vegedaz/) (*read_data_from_Vegedaz_beta_dez_24.r*) *only Beta from Dez23*
    - [Ecological Indicator Values for Europe (EIVE) 1.0](https://vcs.pensoft.net/article/98324/) paper by Dengler et al 2023. (*dengler <- read.xlsx("https://vcs.pensoft.net/article/98324/download/suppl/38/", sheet = 2)*(*require(openxlsx)*))
    - [Veg, a software to manage vegetation data](https://www.maerki.com/maerki_informatik/veg/index.html)  (read_data_from_Veg_2015.r)
      - based on Flora Helvetica v.5 und Flora Indicativa	ca. 2014
- ... transpose your data frame from one survey per column to one survey per row. This is often a prerequisite to do the analysis (*transpose_data.r*)
- ... align your species names with the collected data (*Correct_species_names.r*)
- ... calculate the weighted or mean indicator values for your vegetation surveys (*Calculate_indicator_values.r*)

- There is also an example that illustrates how it should work. (*example.r*)

## Download all scripts

use this link: [https://github.com/bazz-the-spazz/vegebaz/zipball/master/](https://github.com/bazz-the-spazz/vegebaz/zipball/master/)
