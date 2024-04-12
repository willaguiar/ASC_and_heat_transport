# Draft figure list

Note that there is detailed discussion of these figures in [this issue](https://github.com/willaguiar/ASC_and_heat_transport/issues/26). Let's keep the discussion on that issue, and we can use this page as a summary of that discussion.

### NOTE: Each figure should have its own notebook.

## Figure 1 (Paul)

Figure 1 shows where the regimes are, and what the ASC and heat transport look like on average in each regime.

a) Regime map.

<img width="700" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/20695740/6b31501b-1666-4042-b9d4-ffc523ccfc39">

b-c) Average cross slope heat transport and ASC at each depth in each of the regimes.

<img width="600" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/20695740/eb3bf978-e19c-49bb-8c56-884f8c023aa6">

## Figure 2 (Ellie, Taimoor)

Figure 2 shows monthly correlations in each regime. The top panel shows correlations with depth-average ASC and the middle panel shows correlations with the ASC at the local depth. Only show correlations (red lines), not regressions (green lines). Use fixed regime mask and longitude binning (width yet to be decided), with weighted averaging of the longitude bins to account for how many grid points in each bin. On the right of the two top rows, we want to show scatter pdfs for each regime in different colours, only using points from the deeper regions (>800m?). Maybe indicate significance of correlations with bolder lines?

Main points for this figure:
1. There are only weak correlations between the depth average ASC and the cross slope heat transport, with the exception of heat transport at 600-800m depth in the Surface regime. 
2. Correlations are highest at the surface (but heat transport here is low, so not important) and at depth.
3. Difference / similarity between regimes?

Top row: correlations with depth average ASC, middle row: correlations with ASC at each depth (both need pdf scatterplot on right). Bottom row: spatial map of correlation, averaged over deep water column (400-1200m?) at each longitude around Antarctica (but not stereographic projection): Using 5deg bins?
<img width="800" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/54887782/459d9bc3-d7fb-4ec2-9ba6-2efce034716e">

## Figure 3 (Wilton,Taimoor)

Figure 3 shows how correlations vary with different time averaging (interannual variability, seasonal climatology and high frequency (daily minus seasonal climatology)). As for figure 2, only show correlations (red lines), not regressions (green lines). Use fixed regime mask and longitude binning (width yet to be decided), with weighted averaging of the longitude bins to account for how many grid points in each bin. Instead of the scatterplots currently shown on the right, we want pdfs for each regime in different colours, only using points from the deeper regions (>800m?).

Main points here:
1. High frequency correlations (eddy timescales) are surprisingly low.
2. Strong correlations on seasonal and interannual time scales.

<img width="1100" alt="Screenshot 2024-03-28 at 11 00 27 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/1295b13b-6d30-4b5f-8846-4bea4d8f52d1">

## Figure 4 (Adele to try)

Schematic, or transects with real model output, highlighting how the ASC and heat transport directions/magnitude are changing in each regime, just for the depths/regimes where r^2>0.5. e.g.:

<img width="400" alt="Screenshot 2024-01-18 at 10 51 45 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/c2d6fa76-4c32-4c0c-9be4-b399d81290fe">

<img width="600" alt="Screenshot 2024-04-07 at 11 37 59 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/c0b7b911-6f77-4ffb-be57-9b7a7083a861">

## Supplementary figures

- Sensitivity to bin size (repeat of Fig 2 with different bin sizes). **Fabio**
  
- Sensitivity to including zonal and vertical convergence terms. Two panels: a) Repeat of Figure 2b without zonal convergence included, b) Repeat of Figure 2b also with added vertical convergence. **Wilton**

- Time series for sample regions in each regime showing interannual variability of deep ASC and deep heat transport.

