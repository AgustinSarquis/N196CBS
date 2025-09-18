# REQUESTING ALEURITES MOLUCCANUS (KUKUI) TRAITS FROM TRY DATABASE
# load TRYtraits table
traitlist=c(
  "Plant biomass and allometry: Coarse root dry mass per plant",
  "Plant biomass and allometry: Fine root dry mass per plant",
  "Plant biomass and allometry: Leaf dry mass per plant",
  "Plant biomass and allometry: Root dry mass per plant",
  "Wood (heartwood) carbon (C) content per heartwood dry mass",
  "Wood (sapwood) carbon (C) content per dry mass",
  "Twig carbon content per twig dry mass",
  "Stem sapwood carbon (C) content: empirical form factor for calculation"	,
  "Stem carbon (C) content per stem dry mass"	,
  "Root growth: root carbon (C) production per ground area and time",
  "Root exudation: carbon (C) exudation as fraction of root carbon (C) growth",
  "Root exudation: carbon (C) exudation as fraction of root dry mass",
  "Root decomposition: fraction of carbon (C) remaining",
  "Root carbon (C) content per root dry mass",
  "Root carbon (C) content per ground area"	,
  "Plant carbon (C) allocation to root, stem, leaves",
  "Plant carbon (C) allocation: aboveground/belowground net primary production (NPP) ratio",
  "Plant carbon (C) content per plant dry mass",
  "Litter decomposition: fraction of carbon (C) remaining"	,
  "Litter (other than leaf: twigs, bark, flowers…) carbon (C) content per litter dry mass",
  "Litter (leaf) carbon (C) content per leaf litter dry mass"	,
  "Litter (fine root) carbon (C) content per fine root litter dry mass"	,
  "Fine root growth: fine root carbon (C) production per ground area and time"	,
  "Fine root exudation: carbon (C) exudation as fraction of fine root dry mass"	,
  "Fine root decomposition: fraction carbon (C) remaining"	)
trait_ids <-TRYtraits$TraitID[match(traitlist,TRYtraits$Trait)]
