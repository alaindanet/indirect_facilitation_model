## Cool
load(file = "./inst/scenar_test_u=0.Rdata")

double_states<- c("coexistence", "nurse", "protegee", "desert",
  "protegee_desert", "nurse_desert", "coexistence_desert", "protegee_nurse",
  "unkown")
my_colours <- c(coexistence = "orange", nurse = "green", protegee = "black",
  desert = "#C19A6B", protegee_desert = "gray60", 
  nurse_desert = "green2", coexistence_desert = "orange2", protegee_nurse =
    "darkgreen", unkown = "gray")

plot_diagram(u0,
  param = c(x = "b", y = "g"),
  possible_states = double_states,
  col_states = my_colours,
  debug_mode = FALSE, type = "double")
#TODO: make a plot of multiple states
