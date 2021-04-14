library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# Generate nodes for MBNMA classes
classnodes <-
  create_node_df(
    n = 6,
    label = c("data.frame\\l",
              "mb.network\\l - plot()\\l - summary()\\l",
              "mbnma\\l - plot()\\l - summary()\\l",
              "mb.predict\\l - plot()\\l - summary()\\l",
              "nodesplit\\l - plot()\\l - summary()\\l",
              "mb.rank\\l - plot()\\l - summary()\\l"),
    color = "black",
    fontname="Consolas",
    shape = "rectangle",
    fillcolor = "Honeydew",
    fontcolor = "black",
    fontalign="left",
    height=0.7,
    width=2
  )

# Generate edges between classes
# classedges <-
#   create_edge_df(
#     from = c("1", "2", "2", "2", "3", "3"),
#     to =   c("2", "3", "5", "6", "4", "7"),
#     color="black",
#     rel = "a",
#     fontname="Consolas"
#     )

g <- create_graph(nodes_df = classnodes, attr_theme = "tb")

# Generate function nodes
funnodes <-
  create_node_df(
    n=8,
    label=c("mb.network()", "mb.run()", "mb.nodesplit()", "predict()", "rank()",
            "timeplot()", "fitplot()", "devplot()"), #"cumrank()"),
    shape="rectangle",
    fontname="Consolas",
    fillcolor = "white",
    fontcolor = "black",
    color="black",
    width=1.5,
    style="dashed"
  )
g <- add_node_df(g, funnodes)


# Generate edges between classes and functions
funedges <-
  create_edge_df(
    from = c("1", "7", "2", "8", "3", "10", "2", "9", "3", "11", "2", "2", "2"),
    to =   c("7", "2", "8", "3", "10", "4", "9", "5", "11", "6", "12", "13", "14"),
    color="black",
    rel = "a",
    fontname="Consolas",
    arrowhead=c(rep(c("None", "normal"),6), rep("none",4))
  )
g <- add_edge_df(g, funedges)


# Set graph attributes
g <- g %>%
  add_global_graph_attrs( attr = "splines",
                          value = "ortho",
                          attr_type = "graph")

# Render graph
render_graph(g)

# Save graph
render_graph(g) %>% export_svg %>% charToRaw %>%
  rsvg_png("~/MBNMA/MBNMA R Package/Time/MBNMAtime/man/figures/functionstructure.png")



