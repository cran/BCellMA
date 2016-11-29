#' The frequencies of nucleotide mutations and the frequencies of transition and transversion
#'
#' Calculate nucleotide (A, T, C and G) mutations, transition, transversion and their frequency.
#' @param data A colums from 11 to 22 in the IMGT Table 8 caled "8_V-REGION-nt-mutation-statistics.txt".
#' @return Output is four different values:
#' @return vregion_nt_mut number of nucleotide mutations in V region.
#' @return transition_rel the frequencies of transition.
#' @return transverion_rel the frequencies of transversion.
#' @return trans_transv_anzahl number of transition and transversion in V region.
#' @examples
#' data(IMGTtab8)
#' bm_proband<- nucleotide_mutation(IMGTtab8[,11:22])
#' percentlabels<- round(bm_proband$vregion_nt_mut/sum(bm_proband$vregion_nt_mut)*100)
#' pielabels<- paste(percentlabels, "%", sep="")

#' pie(bm_proband$vregion_nt_mut, col=c("grey50","black","grey90","white"),
#'     labels=pielabels, cex=1.5, radius = 0.3)
#' legend("right", c(" prod. IGHV"), cex=1.5)
#' legend("topleft", c("A","G","T","C"), cex=2, fill=c("grey50","black",
#'                                                    "grey90","white"))

#' @export
#'
nucleotide_mutation <- function(data){

  # data      - Daten aus der IMGT-Tabelle 8
  # Ergebnis  - Anzahl an Nukleotidmutationen, Transitionen, Transversionen
  #             und relative Transitionen/Transversionen

  names(data)<-c("V_REGION_a_g", "V_REGION_g_a", "V_REGION_c_t",
                 "V_REGION_t_c", "V_REGION_a_c", "V_REGION_c_a",
                 "V_REGION_a_t", "V_REGION_t_a", "V_REGION_g_c",
                 "V_REGION_c_g", "V_REGION_g_t", "V_REGION_t_g")


  ag<-substr(data$V_REGION_a_g, 1,1)
  sum_ag<-sum(as.numeric(ag))

  tg<-substr(data$V_REGION_t_g, 1,1)
  sum_tg<-sum(as.numeric(tg))

  cg<-substr(data$V_REGION_c_g, 1,1)
  sum_cg<-sum(as.numeric(cg))

  ga<-substr(data$V_REGION_g_a, 1,1)
  sum_ga<-sum(as.numeric(ga))

  ta<-substr(data$V_REGION_t_a, 1,1)
  sum_ta<-sum(as.numeric(ta))

  ca<-substr(data$V_REGION_c_a, 1,1)
  sum_ca<-sum(as.numeric(ca))

  gt<-substr(data$V_REGION_g_t, 1,1)
  sum_gt<-sum(as.numeric(gt))

  ct<-substr(data$V_REGION_c_t, 1,1)
  sum_ct<-sum(as.numeric(ct))

  at<-substr(data$V_REGION_a_t, 1,1)
  sum_at<-sum(as.numeric(at))

  gc<-substr(data$V_REGION_g_c, 1,1)
  sum_gc<-sum(as.numeric(gc))

  ac<-substr(data$V_REGION_a_c, 1,1)
  sum_ac<-sum(as.numeric(ac))

  tc<-substr(data$V_REGION_t_c, 1,1)
  sum_tc<-sum(as.numeric(tc))

  a<-sum_ga+sum_ca+sum_ta
  g<-sum_ag+sum_tg+sum_cg
  t<-sum_ct+sum_gt+sum_at
  c<-sum_ac+sum_gc+sum_tc

  vregion_nt_mut<-c(a,g,t,c)

  transition1<-sum(c(sum_ag, sum_ga, sum_ct, sum_tc))
  transverion1<-sum(c(sum_ac, sum_at, sum_ca, sum_cg,
                      sum_gc, sum_gt, sum_ta, sum_tg))

  trans_transv<-c(transition1, transverion1)

  transition<-(transition1/sum(transition1,transverion1)*100)
  transverion<-(transverion1/sum(transition1,transverion1)*100)

  ergebnis<-list(vregion_nt_mut, transition, transverion, trans_transv)
  names(ergebnis)<-c("vregion_nt_mut", "transition_rel",
                     "transverion_rel", "trans_transv_anzahl")
  return(ergebnis)
}

