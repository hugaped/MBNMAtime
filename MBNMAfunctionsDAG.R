library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# Generate nodes for MBNMA classes
classnodes <-
  create_node_df(
    n = 8,
    label = c("data.frame\\l",
              "mb.network\\l - plot()\\l - summary()\\l",
              "mbnma\\l - plot()\\l - summary()\\l",
              "mb.predict\\l - plot()\\l - summary()\\l",
              "nodesplit\\l - plot()\\l - summary()\\l",
              "mb.rank\\l - plot()\\l - summary()\\l",
              "relative.array\\l - plot()",
              "nma"
              ),
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
    n=12,
    label=c("mb.network()", "mb.run()", "mb.nodesplit()", "predict()", "rank()",
            "timeplot()", "fitplot()", "devplot()", "get.relative()", "cumrank()",
            "nma.run()", "binplot()"),
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
    from = c("1", "9", "2", "10", "3", "12", "2", "11", "3", "13", "2", "2", "2", "3", "17", "6",
             "2", "2", "19", "4"),
    to =   c("9", "2", "10", "3", "12", "4", "11", "5", "13", "6", "14", "15", "16", "17", "7", "18",
             "20", "19", "8", "13"),
    color="black",
    rel = "a",
    fontname="Consolas",
    arrowhead=c("none", "normal", "none", "normal", "none", "normal", "none", "normal",
                "none", "normal", "normal", "normal", "normal", "none", "normal", "normal",
                "normal", "none", "normal", "none")
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



