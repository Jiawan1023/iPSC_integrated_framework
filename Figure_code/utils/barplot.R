library(ggplot2)
library(dplyr)
#input df and group prder
# Define custom order
group_order 

# Convert group to factor with specified order
df$group <- factor(df$group, levels = group_order)

# Recalculate summary with ordered factor
df_summary <- df %>%
  group_by(group) %>%
  summarise(
    mean_value = mean(Nkx2.1, na.rm = TRUE),
    se = sd(Nkx2.1, na.rm = TRUE) / sqrt(n())
  )

# Use same factor levels in summary
df_summary$group <- factor(df_summary$group, levels = group_order)

#input group colors
# Plot
p <- ggplot() +
  geom_bar(
    data = df_summary,
    aes(x = group, y = mean_value),
    stat = "identity",
    color = "black",
    fill = NA,
    width = 0.5
  ) +
  geom_errorbar(
    data = df_summary,
    aes(x = group, ymin = mean_value - se, ymax = mean_value + se),
    width = 0.2
  ) +
  geom_jitter(
    data = df,
    aes(x = group, y = Nkx2.1, color = group),
    width = 0.2,
    size = 1.5,
    alpha = 0.8
  ) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Nkx2.1 Expression by Group",
    x = "Group",
    y = "Nkx2.1 Expression",
    color = "Group"
  )

print(p)


#---------------------------------
library(ggplot2)
library(dplyr)
df <- EdU_summary_input
# Custom group order
#input group_order

# Set factor levels
df$group <- factor(df$group, levels = group_order)
df$chase_hours <- factor(df$chase_hours, levels = c("nochase", "24hrchase", "48hrchase"))

# Custom colors for dots
#input group_colors

# Summarize mean and sd per group and chase_hours
df_summary <- df %>%
  group_by(group, chase_hours) %>%
  summarise(
    mean_value = mean(positive_cells, na.rm = TRUE),
    sd = sd(positive_cells, na.rm = TRUE),
    .groups = "drop"
  )

# Ensure the summary has the same factor levels
df_summary$group <- factor(df_summary$group, levels = group_order)
df_summary$chase_hours <- factor(df_summary$chase_hours, levels = c("nochase", "24hrchase", "48hrchase"))

# Plot
p <- ggplot() +
  geom_bar(
    data = df_summary,
    aes(x = group, y = mean_value, group = chase_hours),
    stat = "identity",
    position = position_dodge(width = 0.8),
    fill = NA,
    color = "black",
    width = 0.6
  ) +
  geom_errorbar(
    data = df_summary,
    aes(x = group, ymin = mean_value - sd, ymax = mean_value + sd, group = chase_hours),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_jitter(
    data = df,
    aes(x = group, y = positive_cells, color = group),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1,
    alpha = 0.8
  ) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~chase_hours) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Positive Cells by Group and Chase Condition",
    x = "Group",
    y = "Positive Cells",
    color = "Group"
  )

print(p)


#width10 height7
