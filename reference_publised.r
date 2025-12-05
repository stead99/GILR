library(ggplot2)
library(data.table)
library(ggpubr)
library(wesanderson)
library(stringr)

### number ---
dt_aff <- fread('/Users/stead/JS_data/wang-lab/report/20251124/supp/Most Relevant Affiliations/MostRelAffiliationsTable.csv')
dt_aff[, Articles := as.numeric(Articles)][, Affiliation := str_to_title(Affiliation)]
pdf(file = '/Users/stead/JS_data/wang-lab/report/20251124/supp/Affiliations.pdf', 6, 3)
pt <- ggdotchart(dt_aff[1:11][Affiliation != 'Electronic Address: '], x = "Affiliation", y = "Articles",
           color = "Red",                                # Color by groups
           sorting = "descending",                        # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,     
           dot.size = 6,        
           label = dt_aff[1:10][['Articles']],
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),
           ggtheme = theme_pubr())
print(pt)
dev.off()

### timing ---
dt_time <- fread('/Users/stead/JS_data/wang-lab/report/20251124/supp/Countries\' Production over Time/CountryOverTimeTable.csv')
dt_time[, Articles := as.numeric(Articles)]
pdf(file = '/Users/stead/JS_data/wang-lab/report/20251124/supp/country_per_year.pdf', 6, 3)
pt <- ggplot(dt_time, aes(x = Year, y = Articles, group = Country)) +
  geom_line(aes(color = Country), size = 1) +
  scale_color_manual(values=c('lightskyblue1', '#EF3B2C', '#41AB5D', '#FED976', '#034E7B')) +
  theme_bw()
print(pt)
dev.off()
