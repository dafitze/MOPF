# Custom ggplot theme to make pretty plots
# Get the font at https://fonts.google.com/specimen/Barlow+Semi+Condensed
theme_clean <- function() {
  theme_minimal() + #base_family = "Barlow Semi Condensed") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90",
                                          linetype = 2),
          text = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(), #family = "BarlowSemiCondensed-Bold"),
          axis.title = element_text(size = 18), #family = "BarlowSemiCondensed-Medium"),
          strip.text = element_text( #family = "BarlowSemiCondensed-Bold",
                                    size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA))
}