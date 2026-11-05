module TaxodistPlotsExt

# Loaded automatically when both Taxodist and Plots are loaded. The plotting
# functions of Taxodist (`plot_taxodist_cluster`, `plot_taxodist_ord`, and
# `taxo_heatmap`) obtain the Plots module from this extension, so that
# `using Taxodist` alone does not load Plots.

using Taxodist
using Plots

end
