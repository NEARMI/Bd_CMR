if (plot_labs == "Spec-Pop") {

spec_pop_plot_labels <- c(
"Ambystoma cingulatum
SMNWR East"

, "Ambystoma cingulatum
SMNWR West"

, "Anaxyrus boreas
Blackrock (Complex)"

, "Anaxyrus boreas
Blackrock (Heron)"

, "Anaxyrus boreas
Jones Pond"

, "Anaxyrus boreas
Sonoma Mountain"
  
, "Anaxyrus boreas
Two Medicine"
  
, "Pseudacris maculata
Lily Pond"
  
, "Pseudacris maculata
Matthews Pond"
  
, "Notophthalmus viridescens
Mud Lake"
  
, "Notophthalmus viridescens
Scotia Barrens"
  
, "Notophthalmus viridescens
SMNWR West"
  
, "Notophthalmus viridescens
Springfield"
  
, "Rana pretiosa
Dilman Meadows"
  
, "Rana boylii
Fox Creek"
  
, "Rana luteiventris
Jones Pond"
  
, "Rana luteiventris
Lost Horse"
  
, "Rana draytonii
San Francisquito"
  
, "Rana sierrae
Summit Meadow"
  
, "Rana cascadae
Three Creeks"
    )

} else if (plot_labs == "Pop") {

if (!fit_ind_mehg) {
  
spec_pop_plot_labels <- c(
  "FL - SMNWR East"
, "FL - SMNWR West"
, "WY - Blackrock (Complex)"
, "WY - Blackrock (Heron)"
, "MT - Jones Pond"
, "CA - Sonoma Mountain"
, "MT - Two Medicine"
, "CO - Lily Pond"
, "CO - Matthews Pond"
, "WI - Mud Lake"
, "PA - Scotia Barrens"
, "FL - SMNWR West"
, "MA - Springfield"
, "OR - Dilman Meadows"
, "CA - Fox Creek"
, "MT - Jones Pond"
, "MT - Lost Horse"
, "CA - San Francisquito"
, "CA - Summit Meadow"
, "OR - Three Creeks"
    )

} else {
  
spec_pop_plot_labels <- c(
  "WY - Blackrock (Heron)"
, "MT - Jones Pond"
, "CA - Sonoma Mountain"
, "CO - Lily Pond"
, "CO - Matthews Pond"
, "OR - Dilman Meadows"
, "CA - Fox Creek"
, "MT - Jones Pond"
, "MT - Lost Horse"
, "CA - San Francisquito"
, "OR - Three Creeks"
    )
  
}

} else {
  print("Not supported"); break
}

if (!fit_ind_mehg) {

spec_labs <- c(
   expression(italic("Ambystoma cingulatum")~"[Intercept]")
,  expression(italic("Anaxyrus boreas"))
,  expression(italic("Pseudacris maculata"))
,  expression(italic("Notophthalmus viridescens"))
,  expression(italic("Rana")~"spp.")
)

spec_names <- c(
  "Ambystoma cingulatum"
, "Anaxyrus boreas"
, "Pseudacris maculata"
, "Notophthalmus viridescens"
, "Rana spp." 
)

} else {
  
spec_labs <- c(
   expression(italic("Anaxyrus boreas"))
,  expression(italic("Pseudacris maculata"))
,  expression(italic("Rana")~"spp.")
)

spec_names <- c(
  "Anaxyrus boreas"
, "Pseudacris maculata"
, "Rana spp." 
) 
  
}

if (!fit_ind_mehg) {

facet_names <- list(
  
  expression(
    atop(
      bolditalic('Ambystoma cingulatum')
      , bold("SMNWR East")
    )
    )
  
, expression(
    atop(
      bolditalic('Ambystoma cingulatum')
    , bold("SMNWR West")
    )
  )

, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Blackrock (Complex)")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Blackrock (Heron)")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Jones Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Sonoma Mountain")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Two Medicine")
    )
  )
  
, expression(
    atop(
      bolditalic('Pseudacris maculata')
    , bold("Lily Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Pseudacris maculata')
    , bold("Matthews Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Notophthalmus viridescens')
    , bold("Mud Lake")
    )
  )
  
, expression(
    atop(
      bolditalic('Notophthalmus viridescens')
    , bold("Scotia Barrens")
    )
  )
  
, expression(
    atop(
      bolditalic('Notophthalmus viridescens')
    , bold("SMNWR West")
    )
  )
  
, expression(
    atop(
      bolditalic('Notophthalmus viridescens')
    , bold("Springfield")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana pretiosa')
    , bold("Dilman Meadows")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana boylii')
    , bold("Fox Creek")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana luteiventris')
    , bold("Jones Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana luteiventris')
    , bold("Lost Horse")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana draytonii')
    , bold("San Francisquito")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana sierrae')
    , bold("Summit Meadow")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana cascadae')
    , bold("Three Creeks")
    )
  )
  
)

} else {
  
facet_names <- list(

 expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Blackrock (Heron)")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Jones Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Anaxyrus boreas')
    , bold("Sonoma Mountain")
    )
  )
  
, expression(
    atop(
      bolditalic('Pseudacris maculata')
    , bold("Lily Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Pseudacris maculata')
    , bold("Matthews Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana pretiosa')
    , bold("Dilman Meadows")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana boylii')
    , bold("Fox Creek")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana luteiventris')
    , bold("Jones Pond")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana luteiventris')
    , bold("Lost Horse")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana draytonii')
    , bold("San Francisquito")
    )
  )
  
, expression(
    atop(
      bolditalic('Rana cascadae')
    , bold("Three Creeks")
    )
  )
  
)
  
}

facet_labeller <- function(variable,value){
  return(facet_names[value])
}

