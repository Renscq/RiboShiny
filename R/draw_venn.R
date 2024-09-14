########################################
# draw_venn.R
# 
# venn_data: A list of vectors, each vector represents a set
# show_intersect: Whether to show the intersection of sets
# set_color: The color of the set
# set_size: The size of the set
# label: The label of the set
# label_alpha: The alpha of the label
# label_geom: The geom of the label
# label_color: The color of the label
# label_size: The size of the label
# label_txtWidth: The width of the label
# edge_lty: The line type of the edge
# edge_size: The size of the edge
# force_upset: Whether to force the upset plot
# nintersects: The number of intersects
# relative_height: The relative height of the plot
# relative_width: The relative width of the plot
# title: The title of the plot


require(ggVennDiagram)
require(UpSetR)

draw_venn <- function(
    venn_data = NULL, 
    show_intersect = FALSE,
    set_color = "black",
    set_size = 5,
    font_size = 12,
    
    label = "both",
    label_alpha = 0.5,
    label_geom = "text",
    label_color = "black",
    label_size = 5,
    label_txtWidth = 40,
    
    edge_lty = "solid",
    edge_size = 1,
    
    low_color = 'white',
    high_color = 'red',
    
    force_upset = FALSE,
    nintersects = 40,
    relative_height = 3,
    relative_width = 0.3
    ) {
  
  # browser()
  
  # Draw Venn Diagram
  venn_plot <- ggVennDiagram(
    venn_data,
    # category.names = names(venn_data),
    show_intersect = show_intersect,
    set_color = set_color,
    set_size = set_size,
    
    label = label,
    label_alpha = label_alpha,
    label_geom = label_geom,
    label_color = label_color,
    label_size = label_size,
    label_percent_digit = 0,
    label_txtWidth = label_txtWidth,
    
    edge_lty = edge_lty,
    edge_size = edge_size,

    force_upset = force_upset,
    nintersects = nintersects,
    order.intersect.by = c("size", "name", "none"),
    order.set.by = c("size", "name", "none"),
    relative_height = relative_height,
    relative_width = relative_width
  ) +
    scale_fill_gradient(low = low_color, high = high_color)

  return(venn_plot)
}
