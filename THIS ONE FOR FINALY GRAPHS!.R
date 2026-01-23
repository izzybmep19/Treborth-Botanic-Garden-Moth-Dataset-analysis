library(mgcv)
library(DBI)
library(dbplyr)
library(dplyr)
library(odbc)
library(rbms)
library(readr)
library(readxl)
library(tidyverse)
library(RODBC)
library(cli)
library(MASS)
library(vegan)
library(emmeans)
library(patchwork)


#######Get data from access into R###########
con <- odbcConnectAccess2007(
  "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Izzys Database.accdb"
)

sp_data = sqlFetch(con, "IB_ALL_DATA_AND_WIND")
View(sp_data)
#######Biodiversity indices and graphs#######
moth_summary <- sp_data %>%
  group_by(Date, `Common Name`) %>%
  summarise(Abundance = sum(COUNT, na.rm = TRUE)) %>%
  ungroup()

moth_summary <- moth_summary %>%
  mutate(Year = lubridate::year(as.Date(Date, format = "%d/%m/%Y")))


moth_yearly <- moth_summary %>%
  group_by(Year, `Common Name`) %>%
  summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

moth_matrix_year <- moth_yearly %>%
  pivot_wider(
    names_from = `Common Name`,
    values_from = Total,
    values_fill = 0
  )

moth_matrix_year$Year <- as.integer(moth_matrix_year$Year)

diversity_df <- moth_matrix_year %>%
  mutate(
    Shannon  = diversity(dplyr::select(., -Year), index = "shannon"),
    Simpson  = diversity(dplyr::select(., -Year), index = "simpson"),
    Richness = specnumber(dplyr::select(., -Year))
  ) %>%
  dplyr::select(Year, Shannon, Simpson, Richness)

ggplot(diversity_df, aes(x = as.numeric(Year), y = Simpson)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Simpson Diversity Index by Year",
       x = "Year",
       y = "Simpson Diversity")+
  xlim(1991,2025)+
  ylim(0,1)

ggplot(diversity_df, aes(x = as.numeric(Year), y = Shannon)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Year",
       x = "Year",
       y = "Shannon Diversity") +
  xlim(1991,2025)+
  ylim(0,5)

diversity_df <- diversity_df %>%
  mutate(Year = factor(Year))


emm_shannon <- emmeans(gam_shannon, ~ Year)

pairs(emm_shannon, adjust = "tukey")

emm_simpson <- emmeans(gam_simpson, ~ Year)

pairs(emm_simpson, adjust = "tukey")

shannon_results <- pairs(emm_shannon, adjust = "tukey")

write.table(
  as.data.frame(shannon_results),
  "shannon_year_comparisons.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
simpson_results <- pairs(emm_simpson, adjust = "tukey")

write.table(
  as.data.frame(simpson_results),
  "simpson_year_comparisons.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#species abundace and richness

Species_Abundance_by_week <- read_excel("C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Access stuff/Results/Species Abundance by week.xlsx")
View(Species_Abundance_by_week)

ggplot(
  Species_Abundance_by_week,
  aes(x = Week, y =`Number Of species`)
) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6),
    method.args = list(family = nb()),
    se = TRUE
  ) +
  facet_wrap(~ Year) +
  theme_minimal()

ggplot(
  Species_Abundance_by_week,
  aes(x = Year, y =`Number Of species`)
) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6),
    method.args = list(family = nb()),
    se = TRUE
  ) +
  theme_minimal()

ggplot(
  Species_Abundance_by_week,
  aes(x = , y =`Mean Species Abundance`)
) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6),
    method.args = list(family = nb()),
    se = TRUE
  ) +
  facet_wrap(~ Year) +
  theme_minimal()


Species_Abundance_by_week <- Species_Abundance_by_week %>%
  rename(
    Number_of_species = `Number Of species`,
    Mean_species_abundance = `Mean Species Abundance`
  )

sink("Abundance_model_sata.txt")

cat("---Species richness---")

gam_number_species <- Species_Abundance_by_week %>%
  group_by(Year) %>%
  do(
    model = gam(
      Number_of_species ~ s(Week, k = 6),
      family = nb(),
      data = .
    )
  )

for (i in seq_along(gam_number_species$model)) {
  cat("Year:", gam_number_species$Year[i], "\n")
  print(summary(gam_number_species$model[[i]]))
  cat("-------------------------------------")
}

# Test overall differences between years
cat("---Year effect on species richness---")

gam_richness_year <- gam(
  Number_of_species ~ Year + s(Week, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

print(summary(gam_richness_year))

#Test whether seasonal patterns differ between years
cat("---Seasonal pattern differences (Richness) ---")

gam_richness_no_interaction <- gam(
  Number_of_species ~ Year + s(Week, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

gam_richness_interaction <- gam(
  Number_of_species ~ Year + s(Week, by = Year, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

print(anova(gam_richness_no_interaction,
            gam_richness_interaction,
            test = "Chisq"))



cat("---Species Abundance---")
gam_mean_abundance <- Species_Abundance_by_week %>%
  group_by(Year) %>%
  do(
    model = gam(
      Mean_species_abundance ~ s(Week, k = 6),
      family = nb(),
      data = .
    )
  )

for (i in seq_along(gam_mean_abundance$model)) {
  cat("Year:", gam_mean_abundance$Year[i], "\n")
  print(summary(gam_mean_abundance$model[[i]]))
  cat("-------------------------------------")
}
#Test overall differences between years
cat("---Year effect on mean species abundance---")

gam_abundance_year <- gam(
  Mean_species_abundance ~ Year + s(Week, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

print(summary(gam_abundance_year))

#Test whether seasonal patterns differ between years
cat("---Seasonal pattern differences (Abundance)---")

gam_abundance_no_interaction <- gam(
  Mean_species_abundance ~ Year + s(Week, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

gam_abundance_interaction <- gam(
  Mean_species_abundance ~ Year + s(Week, by = Year, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

print(anova(gam_abundance_no_interaction,
            gam_abundance_interaction,
            test = "Chisq"))
sink()
#########species abundnace and richness yeras in one graph##############
names(Species_Abundance_by_week)[3] <- "num_species"
names(Species_Abundance_by_week)[5] <- "Abun_species"

ggplot(
  Species_Abundance_by_week,
  aes(x = Year, y =`num_species`)
) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6),
    method.args = list(family = nb()),
    se = TRUE
  ) +
  theme_minimal()
ggplot(
  Species_Abundance_by_week,
  aes(x = Year, y =`Abun_species`)
) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6),
    method.args = list(family = nb()),
    se = TRUE
  ) +
  theme_minimal()

gam_mod <- gam(
  `num_species` ~ s(Year, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

summary(gam_mod)

gam_mod2 <- gam(
  `Abun_species` ~ s(Year, k = 6),
  family = nb(),
  data = Species_Abundance_by_week
)

summary(gam_mod2)
##############Specie Richness and Abundance Vs Climate##############
SpAbLine = sqlFetch(con, "SpAbLine")
view(SpAbLine)
Output_Dir3 = "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Attempt nine"

AvTempAbun = ggplot(SpAbLine,aes(x = AvTemp, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

RainAbun = ggplot(SpAbLine,aes(x = Rain, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

AfDaysAbun = ggplot(SpAbLine,aes(x = AfDays, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

MoonAbun = ggplot(SpAbLine,aes(x = Moon, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

#####Wind speed and direction########
wind_df <- SpAbLine %>%
  dplyr::select(Year, Wind_D, Wind_S, Abundance) %>%
  dplyr::filter(
    !is.na(Year),
    !is.na(Wind_D),
    !is.na(Wind_S),
    !is.na(Abundance)
  )


dir_levels <- c(
  "N", "NNE", "NE", "ENE",
  "E", "ESE", "SE", "SSE",
  "S", "SSW", "SW", "WSW",
  "W", "WNW", "NW", "NNW"
)

wind_df <- wind_df %>%
  mutate(
    Wind_D = ifelse(Wind_D == 360, 0, Wind_D),
    wind_dir = cut(
      Wind_D,
      breaks = seq(0, 360, by = 22.5),
      include.lowest = TRUE,
      right = FALSE,
      labels = dir_levels
    )
  ) %>%
  filter(!is.na(wind_dir))


wind_df <- wind_df %>%
  mutate(
    wind_speed_bin = cut(
      Wind_S,
      breaks = c(0, 2, 5, 8, 12, Inf),
      labels = c("0–2", "2–5", "5–8", "8–12", "12+"),
      include.lowest = TRUE
    )
  )

wind_sum_all <- wind_df %>%
  group_by(wind_dir, wind_speed_bin) %>%
  summarise(
    total_Abun = mean(Abundance),
    .groups = "drop"
  )

## Ensure correct direction order
wind_sum_all$wind_dir <- factor(
  wind_sum_all$wind_dir,
  levels = dir_levels
)

## Ensure all directions exist (even zero-count)
wind_sum_all <- wind_sum_all %>%
  tidyr::complete(
    wind_dir,
    wind_speed_bin,
    fill = list(total_Abun = 0)
  )

WindAbun = ggplot(
  wind_sum_all,
  aes(
    x = wind_dir,
    y = total_Abun,
    fill = wind_speed_bin
  )
) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  coord_polar(start = -pi/16) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_fill_viridis_d(name = "Wind speed") +
  labs(
    title = "Wind Rose (All Years)",
    y = "Total abundance"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

#################

WeekAbun = ggplot(SpAbLine,aes(x = Week, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

MonthAbun = ggplot(SpAbLine,aes(x = Month, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

YearAbun = ggplot(SpAbLine,aes(x = Year, y =Abundance)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

AbunCombined_plot <- (WeekAbun | MonthAbun | YearAbun) / 
  (AvTempAbun | MoonAbun | RainAbun) / 
  (AfDaysAbun | WindAbun)+
  plot_annotation(
    title = paste("All Trends for Abundance"),
    tag_levels = "A"
  )
ggsave(
  filename = file.path(Output_Dir3, paste0("All Trends for Abundance.jpeg")),
  plot = AbunCombined_plot,
  width = 15,
  height = 20,
  device = ragg::agg_jpeg
)
##rich##
AvTempRich = ggplot(SpAbLine,aes(x = AvTemp, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth") +
  ylim(0,120)

RainRich = ggplot(SpAbLine,aes(x = Rain, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")
AfDaysRich = ggplot(SpAbLine,aes(x = AfDays, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")+
  ylim(0,150)

MoonRich = ggplot(SpAbLine,aes(x = Moon, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

#####Wind Wspeed and direction########
wind_df <- SpAbLine %>%
  dplyr::select(Year, Wind_D, Wind_S, Richness) %>%
  dplyr::filter(
    !is.na(Year),
    !is.na(Wind_D),
    !is.na(Wind_S),
    !is.na(Richness)
  )


dir_levels <- c(
  "N", "NNE", "NE", "ENE",
  "E", "ESE", "SE", "SSE",
  "S", "SSW", "SW", "WSW",
  "W", "WNW", "NW", "NNW"
)

wind_df <- wind_df %>%
  mutate(
    Wind_D = ifelse(Wind_D == 360, 0, Wind_D),
    wind_dir = cut(
      Wind_D,
      breaks = seq(0, 360, by = 22.5),
      include.lowest = TRUE,
      right = FALSE,
      labels = dir_levels
    )
  ) %>%
  filter(!is.na(wind_dir))


wind_df <- wind_df %>%
  mutate(
    wind_speed_bin = cut(
      Wind_S,
      breaks = c(0, 2, 5, 8, 12, Inf),
      labels = c("0–2", "2–5", "5–8", "8–12", "12+"),
      include.lowest = TRUE
    )
  )

wind_sum_all <- wind_df %>%
  group_by(wind_dir, wind_speed_bin) %>%
  summarise(
    total_Rich = mean(Richness),
    .groups = "drop"
  )

## Ensure correct direction order
wind_sum_all$wind_dir <- factor(
  wind_sum_all$wind_dir,
  levels = dir_levels
)

## Ensure all directions exist (even zero-count)
wind_sum_all <- wind_sum_all %>%
  tidyr::complete(
    wind_dir,
    wind_speed_bin,
    fill = list(total_Rich = 0)
  )

WindRich = ggplot(
  wind_sum_all,
  aes(
    x = wind_dir,
    y = total_Rich,
    fill = wind_speed_bin
  )
) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  coord_polar(start = -pi/16) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_fill_viridis_d(name = "Wind speed") +
  labs(
    title = "Wind Rose (All Years)",
    y = "Total Richness"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
################################

WeekRich = ggplot(SpAbLine,aes(x = Week, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

MonthRich = ggplot(SpAbLine,aes(x = Month, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

YearRich = ggplot(SpAbLine,aes(x = Year, y =Richness)) +
  geom_point(alpha = 0.4, size = 0.8) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    geom = "smooth")

RichCombined_plot <- (WeekRich | MonthRich | YearRich) / 
  (AvTempRich | MoonRich | RainRich) / 
  (AfDaysRich | WindRich)+
  plot_annotation(
    title = paste("All Trends for Richness"),
    tag_levels = "A"
  )
ggsave(
  filename = file.path(Output_Dir3, paste0("All Trends for Richness.jpeg")),
  plot = RichCombined_plot,
  width = 15,
  height = 20,
  device = ragg::agg_jpeg
)
#########Abun vs Climate stats###########
AvTempAbunStat = lm(Abundance~AvTemp, SpAbLine)
RainAbunStat = lm(Abundance~Rain, SpAbLine)
AfDaysAbunStat = lm(Abundance~AfDays, SpAbLine)
MoonAbunStat = lm(Abundance~Moon, SpAbLine)
WindDAbunStat = lm(Abundance~Wind_D, SpAbLine)
WindSAbunStat = lm(Abundance~Wind_S, SpAbLine)
WeekAbunStat = lm(Abundance~Week, SpAbLine)
MonthAbunStat = lm(Abundance~Month, SpAbLine)
YearAbunStat = lm(Abundance~Year, SpAbLine)

AbunModels = list(
  summary(AvTempAbunStat),
  summary(RainAbunStat),
  summary(AfDaysAbunStat),
  summary(MoonAbunStat),
  summary(WindSAbunStat),
  summary(WindDAbunStat),
  summary(WeekAbunStat),
  summary(MonthAbunStat),
  summary(YearAbunStat)
)


sink(file.path(Output_Dir3, "Abundance Models Stats Comprihensive.txt"))

lapply(AbunModels, function(model_list) {
  lapply(model_list, summary)
})

sink()

AbunModels <- list(
  AvTemp  = AvTempAbunStat,
  Rain    = RainAbunStat,
  AfDays  = AfDaysAbunStat,
  Moon    = MoonAbunStat,
  Wind_S  = WindSAbunStat,
  Wind_D  = WindDAbunStat,
  Week    = WeekAbunStat,
  Month   = MonthAbunStat,
  Year    = YearAbunStat
)
sink(file.path(Output_Dir3, "Abundance Models Stats.txt"))

lapply(names(AbunModels), function(name) {
  
  cat("\n====================================\n")
  cat("Predictor:", name, "\n")
  cat("====================================\n\n")
  
  print(summary(AbunModels[[name]]))
  cat("\n\n")
})

sink()
######### Rich vs Climate stats#########
AvTempRichStat = lm(Richness~AvTemp, SpAbLine)
RainRichStat = lm(Richness~Rain, SpAbLine)
AfDaysRichStat = lm(Richness~AfDays, SpAbLine)
MoonRichStat = lm(Richness~Moon, SpAbLine)
WindDRichStat = lm(Richness~Wind_D, SpAbLine)
WindSRichStat = lm(Richness~Wind_S, SpAbLine)
WeekRichStat = lm(Richness~Week, SpAbLine)
MonthRichStat = lm(Richness~Month, SpAbLine)
YearRichStat = lm(Richness~Year, SpAbLine)

RichModels = list(
  summary(AvTempRichStat),
  summary(RainRichStat),
  summary(AfDaysRichStat),
  summary(MoonRichStat),
  summary(WindSRichStat),
  summary(WindDRichStat),
  summary(WeekRichStat),
  summary(MonthRichStat),
  summary(YearRichStat)
)


sink(file.path(Output_Dir3, "Richness Models Stats Comprihensive.txt"))

lapply(RichModels, function(model_list) {
  lapply(model_list, summary)
})

sink()

RichModels <- list(
  AvTemp  = AvTempRichStat,
  Rain    = RainRichStat,
  AfDays  = AfDaysRichStat,
  Moon    = MoonRichStat,
  Wind_S  = WindSRichStat,
  Wind_D  = WindDRichStat,
  Week    = WeekRichStat,
  Month   = MonthRichStat,
  Year    = YearRichStat
)
sink(file.path(Output_Dir3, "Richness Models Stats.txt"))

lapply(names(RichModels), function(name) {
  
  cat("\n====================================\n")
  cat("Predictor:", name, "\n")
  cat("====================================\n\n")
  
  print(summary(RichModels[[name]]))
  cat("\n\n")
})

sink()
#####################Catagory Counts########
SpAbBox = sqlFetch(con, "SpAbBox")
view(SpAbBox)
Output_Dir4 = "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Attempt nine"


AvSize = sqlFetch(con,"AvSize")
BroodNumber = sqlFetch(con,"BroodNumber")
MothCol = sqlFetch(con,"Colour")
Ditribution = sqlFetch(con,"Distribution")
Family = sqlFetch(con,"Family")
FoodCat = sqlFetch(con,"FoodCatagories")
LarvalFoodPlants = sqlFetch(con,"LarvalFoodPlants")
Overwinter = sqlFetch(con,"Overwinter")
Status = sqlFetch(con,"Status")


MothColBar = ggplot(MothCol, aes(x = Colour, y = Count, fill = Colour)) +
  geom_col(color = "grey30", width = 0.8) +
  scale_fill_manual(values = c(
    "black","brown","lightyellow","gold","green","darkgrey",
    "orange","pink","purple","red","grey","white","white",
    "white","yellow"
  )) +
  labs(
    title = "Moth Colours Count",
    x = "Colour",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    text = element_text(size = 30)
    
  )

AvSizeBar = ggplot(AvSize, aes(x = AvSize, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by AvSize",
    x = "AvSize",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

BroodNumBar = ggplot(BroodNumber, aes(x = `Brood Number`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Brood Number",
    x = "Brood Number",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

DistBar = ggplot(Ditribution, aes(x = distribution, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by distribution",
    x = "Distribution",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

FamilyBar = ggplot(Family, aes(x = `Family`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Family",
    x = "Family",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

FoodCatBar = ggplot(FoodCat, aes(x = `Food Catagories`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Food Catagories",
    x = "Food Catagories",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

LarFoodPlantBar = ggplot(LarvalFoodPlants, aes(x = `Larvel Food Plants`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Larval Food Plants",
    x = "Larval Food Plants",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

Overwinterbar = ggplot(Overwinter, aes(x = `Overwinter`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Overwinter",
    x = "Overwinter",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

StatusBar = ggplot(Status, aes(x = `Status`, y = Count)) +
  geom_col() +
  scale_fill_viridis_c() +
  labs(
    title = "Moth Counts by Status",
    x = "Status",
    y = "Count",
    fill = "Average size"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 30))

BarCombi_plot <- (MothColBar | AvSizeBar) / 
  (DistBar | FamilyBar) +
  plot_annotation(
    title = paste("All Barplots For Catagories"),
    tag_levels = "A"
  )
ggsave(
  filename = file.path(Output_Dir4, paste0("All Barplots For Catagories.jpeg")),
  plot = BarCombi_plot,
  width = 30,
  height = 30,
  device = ragg::agg_jpeg
)
BarCombi_plot2 <- (BroodNumBar| Overwinterbar) / 
  ( FoodCatBar|StatusBar) +
  plot_annotation(
    title = paste("All Barplots For Catagories"),
    tag_levels = "A"
  )
ggsave(
  filename = file.path(Output_Dir4, paste0("All Barplots For Catagories2.jpeg")),
  plot = BarCombi_plot2,
  width = 30,
  height = 30,
  device = ragg::agg_jpeg
)
########GLM graphs and stats#######

# Directory into which results of current analysis will saved
output_path = "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Attempt nine"

spp_codes = c(2091, 2189, 2164, 2110, 1773, 1902, 1913, 1681, 1305, 2128, 2318, 1293, 2241, 1727, 2492, 1922, 834, 2127, 1288, 1957, 2237, 970, 2240, 1872, 1766, 2160, 1881, 2489, 2422, 2245, 2089, 2111, 1979, 18, 1646, 1853, 1926, 2484, 2263, 1708, 1702, 2134, 1893, 2118, 1809, 2264, 2030, 2133, 1345, 1883, 2341, 1653, 1762, 2102, 1784, 1761, 1962, 2248, 1887, 1958, 2297, 2155, 1043, 2150, 2443, 2441, 2087, 2176, 1912, 2117, 2232, 1833, 1954, 1945, 2166, 2188, 2187, 17, 2291, 1934, 2243, 1917, 1728, 2444, 1807, 2138, 2254, 1768, 2469, 2289, 2078, 2109, 1338, 874, 2173, 1036, 663, 2063, 1930, 2389, 2028, 1931, 2179, 2139, 2256, 1789, 1920, 2047, 1948, 1376, 2182, 1803, 2350, 2198, 2477, 1769, 1747, 2380, 1750, 1937, 1800, 1799, 2267, 2262, 2069, 2011, 2281, 1903, 1907, 1852, 2061, 2434, 1755, 1651, 2327, 1738, 2050, 1076, 1764, 1524, 1956, 2008, 1693, 2321, 2449, 873, 2130, 1640, 2330, 1390, 1960, 1947, 1923, 1722, 2265, 2353, 1817, 2364, 1720, 2439, 2136, 1951, 1837, 2190, 1777, 2107, 2112, 2322, 1910, 2270, 2280, 1819, 2387, 2299, 1385, 2271, 15, 1652, 1754, 1767, 2278, 1919, 2064, 1398, 2274, 1808, 1940, 1904, 2305, 2450, 1932, 1821, 1936, 1392, 2119, 1921, 1659, 2154, 1802) # replace with your species codes

# ---- safe_nb: only returns a fitted glm.nb object or NULL ----
safe_nb <- function(formula, data, min_rows = 10){
  mf <- try(model.frame(formula, data), silent = TRUE)
  if(inherits(mf, "try-error")) return(NULL)
  df <- mf[complete.cases(mf), ]
  if(nrow(df) < min_rows) return(NULL)
  preds <- all.vars(formula)[-1]
  for(v in preds){
    if(!(v %in% names(df)) || length(unique(df[[v]])) < 2) return(NULL)
  }
  fit <- try(MASS::glm.nb(formula, data = df), silent = TRUE)
  if(inherits(fit, "try-error")) return(NULL)
  if(!inherits(fit, "negbin")) return(NULL)
  return(fit)
}

for(sp in spp_codes){
    df <- subset(sp_data, sp_data$SPECIES == sp)
  n_years <- length(unique(df$YEAR))
  n_weeks <- length(unique(df$WEEKNO))
  
  if (n_years < 4 | n_weeks < 5) {
    cat("Skipping species:", sp,
        "— insufficient data (years:", n_years,
        ", weeks:", n_weeks, ")\n\n")
    next
  }
  
  # --- CREATE PLOTS ---
  p_week <- ggplot(df, aes(x = WEEKNO, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Seasonal (week)") +
    xlim(1, 52)
  
  p_year <- ggplot(df, aes(x = YEAR, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Yearly Trend")+
    xlim(1991,2025)
  
  p_month <- ggplot(df, aes(x = MONTHNO, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Monthly Trend")+
    xlim(1,12)
  
  p_AvTemp <- ggplot(df, aes(x = AvTemp, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Temperature Max")+
    xlim(0,30)
  
  p_Rain <- ggplot(df, aes(x = Rain, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Rain")+
    xlim(0,300)
  
  p_Moon <- ggplot(df, aes(x = Moon, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Moon")+
    xlim(0,1)
  
####wind#####
  wind_df <- df %>%
    dplyr::select(YEAR, WIND_D, WIND_S, COUNT) %>%
    dplyr::filter(
      !is.na(YEAR),
      !is.na(WIND_D),
      !is.na(WIND_S),
      !is.na(COUNT)
    )
  
  
  dir_levels <- c(
    "N", "NNE", "NE", "ENE",
    "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW",
    "W", "WNW", "NW", "NNW"
  )
  
  wind_df <- wind_df %>%
    mutate(
      WIND_D = ifelse(WIND_D == 360, 0, WIND_D),
      wind_dir = cut(
        WIND_D,
        breaks = seq(0, 360, by = 22.5),
        include.lowest = TRUE,
        right = FALSE,
        labels = dir_levels
      )
    ) %>%
    filter(!is.na(wind_dir))
  
  
  wind_df <- wind_df %>%
    mutate(
      wind_speed_bin = cut(
        WIND_S,
        breaks = c(0, 2, 5, 8, 12, Inf),
        labels = c("0–2", "2–5", "5–8", "8–12", "12+"),
        include.lowest = TRUE
      )
    )
  
  wind_sum_all <- wind_df %>%
    group_by(wind_dir, wind_speed_bin) %>%
    summarise(
      total_Rich = sum(COUNT),
      .groups = "drop"
    )
  
  ## Ensure correct direction order
  wind_sum_all$wind_dir <- factor(
    wind_sum_all$wind_dir,
    levels = dir_levels
  )
  
  ## Ensure all directions exist (even zero-count)
  wind_sum_all <- wind_sum_all %>%
    tidyr::complete(
      wind_dir,
      wind_speed_bin,
      fill = list(total_Rich = 0)
    )
  
  p_Wind = ggplot(
    wind_sum_all,
    aes(
      x = wind_dir,
      y = total_Rich,
      fill = wind_speed_bin
    )
  ) +
    geom_col(width = 1, color = "white", linewidth = 0.2) +
    coord_polar(start = -pi/16) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    scale_fill_viridis_d(name = "Wind speed") +
    labs(
      title = "Wind Rose (All Years)",
      y = "Total Count"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
#############
  
  p_Af <- ggplot(df, aes(x = AfDays, y = COUNT)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
                method.args = list(family = nb()), se = TRUE) +
    ggtitle("Air Frost Days")+
    xlim(0,20)
  
  combined_plot <- (p_week | p_month | p_year) / 
    (p_AvTemp | p_Moon | p_Rain) / 
    (p_Af | p_Wind)+
    plot_annotation(
      title = paste("All Trends for", sp),
      tag_levels = "A"
    )
  
  # Save the combined plot
  ggsave(
    file.path(output_path, paste0("combined_plot_", sp, ".png")),
    combined_plot,
    width = 15,
    height = 20
  )
  
  cat("\n--- Finished plotting for species:", sp, "---\n\n")
  
  models <- list(
    Week        = safe_nb(COUNT ~ WEEKNO, df),
    Month       = safe_nb(COUNT ~ MONTHNO, df),
    Year        = safe_nb(COUNT ~ YEAR, df),
    AVTEMP     = safe_nb(COUNT ~ AvTemp, df),
    RAIN        = safe_nb(COUNT ~ Rain, df),
    MOON        = safe_nb(COUNT ~ Moon, df),
    Direction   = safe_nb(COUNT ~ WIND_D, df),
    Speed       = safe_nb(COUNT ~ WIND_S, df),
    AFDAYS      = safe_nb(COUNT ~ AfDays, df),
    
    WeekAv     = safe_nb(COUNT ~ AvTemp + WEEKNO, df),
    WeekRain    = safe_nb(COUNT ~ Rain + WEEKNO, df),
    WeekYear    = safe_nb(COUNT ~ YEAR + WEEKNO, df),
    WeekMonth   = safe_nb(COUNT ~ MONTHNO + WEEKNO, df),
    
    MonthAv    = safe_nb(COUNT ~ MONTHNO + AvTemp, df),
    MonthRain   = safe_nb(COUNT ~ MONTHNO + Rain, df),
    MonthYear   = safe_nb(COUNT ~ MONTHNO + YEAR, df),
    
    YearAv     = safe_nb(COUNT ~ YEAR + AvTemp, df),
    YearRain    = safe_nb(COUNT ~ YEAR + Rain, df),
    
    RainAv     = safe_nb(COUNT ~ Rain + AvTemp, df),
    
    MoonAv     = safe_nb(COUNT ~ Moon + AvTemp, df),
    MoonRain    = safe_nb(COUNT ~ Moon + Rain, df),
    MoonYear    = safe_nb(COUNT ~ Moon + YEAR, df),
    MoonMonth   = safe_nb(COUNT ~ Moon + MONTHNO, df),
    MoonWeek    = safe_nb(COUNT ~ Moon + WEEKNO, df),
    
    DirectionAv  = safe_nb(COUNT ~ WIND_D + AvTemp, df),
    DirectionRain = safe_nb(COUNT ~ WIND_D + Rain, df),
    DirectionYear = safe_nb(COUNT ~ WIND_D + YEAR, df),
    DirectionMonth= safe_nb(COUNT ~ WIND_D + MONTHNO, df),
    DirectionWeek = safe_nb(COUNT ~ WIND_D + WEEKNO, df),
    DirectionMoon = safe_nb(COUNT ~ WIND_D + Moon, df),
    
    SpeedAv   = safe_nb(COUNT ~ WIND_S + AvTemp, df),
    SpeedRain  = safe_nb(COUNT ~ WIND_S + Rain, df),
    SpeedYear  = safe_nb(COUNT ~ WIND_S + YEAR, df),
    SpeedMonth = safe_nb(COUNT ~ WIND_S + MONTHNO, df),
    SpeedWeek  = safe_nb(COUNT ~ WIND_S + WEEKNO, df),
    SpeedMoon  = safe_nb(COUNT ~ WIND_S + Moon, df),
    SpeedDirection = safe_nb(COUNT ~ WIND_S + WIND_D, df),
    
    AfAv   = safe_nb(COUNT ~ AfDays + AvTemp, df),
    AfRain  = safe_nb(COUNT ~ AfDays + Rain, df),
    AfYear  = safe_nb(COUNT ~ AfDays + YEAR, df),
    AfMonth = safe_nb(COUNT ~ AfDays + MONTHNO, df),
    AfWeek  = safe_nb(COUNT ~ AfDays + WEEKNO, df),
    AfMoon  = safe_nb(COUNT ~ AfDays + Moon, df),
    AfDirection = safe_nb(COUNT ~ AfDays + WIND_D, df),
    AfSpeed = safe_nb(COUNT ~ AfDays + WIND_S, df)
  )
  
  # ---- keep only fitted glm.nb objects ----
  valid_models <- Filter(function(x) inherits(x, "negbin"), models)
  
  aic_table <- NULL
  
  if (length(valid_models) > 0) {
    
    aic_vals <- sapply(valid_models, function(m) {
      if (inherits(m, "negbin")) AIC(m) else NA_real_
    })
    
    aic_table <- data.frame(
      Model = names(aic_vals),
      AIC   = as.numeric(aic_vals),
      stringsAsFactors = FALSE
    )
    
    aic_table <- aic_table[order(aic_table$AIC), ]
    
    write.table(aic_table,
                file = file.path(output_path, paste0("AIC_", sp, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  } else {
    message("No valid NB models for species ", sp)
  }
  
  # ---- save summaries of valid models ----
  sink(file.path(output_path, paste0("summary_", sp, "_LM_AIC.txt")))
  cat("=== Species:", sp, "===\n\n")
  if(length(valid_models) > 0){
    cat("AIC table:\n")
    print(aic_table)
    cat("\n\nModel summaries:\n")
    for(nm in names(valid_models)){
      cat("\n--- ", nm, " ---\n")
      print(summary(valid_models[[nm]]))
    }
  } else {
    cat("No valid NB models to summarize.\n")
  }
  sink()
  
} # end species loop


output_dir1 = "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Attempt nine/Overlapped"
output_dir2 = "C:/Users/izzyb/OneDrive/Desktop/uni/Year 4 - MRes/Yes-MSc/Dissertation/Data/Whole time/Attempt nine/Separate"

for(sp in spp_codes){
  Practicemoth = subset(sp_data, sp_data$SPECIES == sp)
  if (nrow(Practicemoth) == 0) next
  
  Practicemoth_filt <- Practicemoth %>%
    filter(YEAR >= 1991, YEAR <= 2025) %>%
    group_by(YEAR) %>%
    filter(
      n() >= 10,                       # enough rows
      n_distinct(WEEKNO) >= 5,         # WEEKNO varies
      var(COUNT, na.rm = TRUE) > 0     # COUNT varies
    ) %>%
    ungroup()
  
  if (nrow(Practicemoth_filt) == 0) next
  
  AllModels <- ggplot(
    Practicemoth_filt,
    aes(x = WEEKNO, y = COUNT, group = YEAR, color = factor(YEAR))
  ) +
    geom_point(alpha = 0.3, size = 0.8) +
    stat_smooth(
      method = "gam",
      formula = y ~ s(x, k = 6),
      method.args = list(family = nb()),
      se = FALSE
    ) +
    guides(color = "none") +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_dir1, paste0(sp, "_GAM_overlayed.jpeg")),
    plot = AllModels,
    width = 10,
    height = 6,
    dpi = 300,
    device = ragg::agg_jpeg
  )
  
  AllModelsSeparate <- ggplot(
    Practicemoth_filt,
    aes(x = WEEKNO, y = COUNT)
  ) +
    geom_point(alpha = 0.4, size = 0.8) +
    stat_smooth(
      method = "gam",
      formula = y ~ s(x, k = 6),
      method.args = list(family = nb()),
      se = TRUE
    ) +
    facet_wrap(~ YEAR) +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_dir2, paste0(sp, "_GAM_Separate.jpeg")),
    plot = AllModelsSeparate,
    width = 14,
    height = 10,
    dpi = 300,
    device = ragg::agg_jpeg
  )
}

