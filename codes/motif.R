enriched.motif <- get.enriched.motif(data = mae,
                                     min.motif.quality = "DS",
                                     probes = unique(Hypo.pair$Probe),
                                     dir.out = "Results_hypo",
                                     label = "hypo",
                                     min.incidence = 10,
                                     lower.OR = 1.1)
