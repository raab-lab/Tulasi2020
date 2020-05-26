# Script of thematic info for the paper

# Population Orders
# Left - Neg, SL, Low, High # 

# Colorscheme

sox9_cols <- c('#636363', '#e5f5f9', '#99d8c9', '#2ca25f')

# theme

theme_sox9 <- function() { 
  theme_minimal() + 
  theme(axis.text = element_text(size = 12, color = 'grey20'),  
        axis.title = element_text(size = 14, color = 'grey10'), 
        strip.text = element_text(size = 16, color = 'grey10'), 
        panel.grid = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"))
} 
