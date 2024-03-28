# Draft figure list

Note that there is detailed discussion of these figures in [this issue](https://github.com/willaguiar/ASC_and_heat_transport/issues/26). Let's keep the discussion on that issue, and we can use this page as a summary of that discussion.

### NOTE: Each figure should have its own notebook.

## Figure 1 (Paul)

Figure 1 shows where the regimes are, and what the ASC and heat transport look like on average in each regime.

a) Regime map.

<img width="700" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/20695740/6b31501b-1666-4042-b9d4-ffc523ccfc39">

b-c) Average cross slope heat transport and ASC at each depth in each of the regimes.

<img width="600" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/20695740/eb3bf978-e19c-49bb-8c56-884f8c023aa6">

## Figure 2 (Ellie)

Figure 2 shows monthly correlations in each regime. The top panel shows correlations with depth-average ASC and the middle panel shows correlations with the ASC at the local depth. Only show correlations (red lines), not regressions (green lines). Use fixed regime mask and longitude binning (width yet to be decided), with weighted averaging of the longitude bins to account for how many grid points in each bin. On the right of the two top rows, we want to show scatter pdfs for each regime in different colours, only using points from the deeper regions (>800m?).

Main points for this figure:
1. There are only weak correlations between the depth average ASC and the cross slope heat transport, with the exception of heat transport at 600-800m depth in the Surface regime. 
2. Correlations are highest at the surface (but heat transport here is low, so not important) and at depth.
3. Difference / similarity between regimes?

Top row depth average ASC (needs pdf scatter on right):
<img width="800" alt="Screenshot 2024-01-18 at 11 10 21 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/88fd7099-f46e-4d2c-83d3-ab776b0780af">

Middle row correlations with ASC at each depth (needs pdf scatter on right):
<img width="800" alt="Screenshot 2024-03-28 at 10 32 11 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/19beff77-b2c6-4690-9df4-27920551ca73">

Bottom row: spatial map of correlation, averaged over deep water column (800-1000?) at each longitude around Antarctica (but not stereographic projection): Using 5deg bins?

<img width="408" alt="Screenshot 2024-03-28 at 10 54 05 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/67ec483f-964d-424a-8749-26c22c28c9c1">

## Figure 3 (Wilton)

Figure 3 shows how correlations vary with different time averaging (interannual variability, seasonal climatology and high frequency (daily minus seasonal climatology)). 

<img width="800" alt="Screenshot 2023-11-23 at 11 00 11 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/3c2e32ef-725b-4f19-961f-1fa0a6fabfdf">

We also want to add some panels showing scatterplots for only the depths where correlations > 0.5. This will help the reader understand what direction the ASC and CSHT are, e.g.

<img width="700" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/70033934/bbcd0bf5-fe33-4082-b2fe-5f40ff1fd2a5">

Maybe we want to include subseasonal correlations here also, using daily data? If we do that, maybe we could move the top row here up to join Figure 2 (so then Figure 2 would show depth avg vs not depth avg, and Figure 3 would just focus on correlations at different time frequencies).

## Figure 4

Figure 4 possibilities:

A) Repeat the above plot for a single basin in each regime, e.g. Totten for Surface, Ross for Deep, Wilkins for Reverse regime.

B) Schematic highlighting how the ASC and heat transport directions/magnitude are changing in each regime, just for the depths/regimes where r^2>0.5. e.g.:

<img width="400" alt="Screenshot 2024-01-18 at 10 51 45 am" src="https://github.com/willaguiar/ASC_and_heat_transport/assets/8506963/c2d6fa76-4c32-4c0c-9be4-b399d81290fe">

Possibly we want both these figures, so 5 total?

## Supplementary figures

- If we use time-varying regimes, show how they vary in time.

- Correlations/regressions in each regime with eddy/mean decomposition.
