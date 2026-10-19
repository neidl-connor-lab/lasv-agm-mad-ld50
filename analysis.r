#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggsurvfit))
options(stringsAsFactors=FALSE)
theme_set(ggpubr::theme_pubr() + 
          theme(legend.position="right",
                legend.background=element_blank()))

# helper variables
cols.outcome <- c("Survived"="black",
                  "Succumbed"="#e41a1c")
cols.group <- c(Succumbed="#e41a1c",
                Symptomatic="black",
                Asymptomatic="darkgrey")
breaks.dates <- c(0, 4, 7, 10, 14, 21, 28, 35)

# load the full Excel table 
sheets <- readxl::excel_sheets("data.xlsx")
df <- lapply(sheets, function(i) { readxl::read_excel("data.xlsx", sheet=i)})
names(df) <- sheets
rm(sheets)

## the basics: survival, scores, and clinical signs ----------------------------
# format outcome and group. order NHPs by day of death and then group
df$animals <- df$animals %>%
              mutate(Outcome=factor(Outcome, levels=names(cols.outcome)),
                     Group=factor(Group, levels=names(cols.group))) %>%
              arrange(Group, Death) %>%
              mutate(NHP=factor(NHP, levels=NHP))

# survival: KM curve 
df$animals %>%
  mutate(Outcome=as.integer(Outcome)-1) %>%
  survfit2(Surv(Death, Outcome) ~ 1, data=.) %>%
  ggsurvfit(theme=ggpubr::theme_pubr(), 
            linewidth=1, lineend="round") +
  scale_y_continuous(limits=c(0, NA), 
                     breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks=breaks.dates) +
  labs(x="Days postinfection",
       y="Survival probability", 
       title="Survival") +
  theme(legend.position="right")
ggsave("analysis/survival.png", units="in", width=4, height=2.25)

# mean and 95% CI TTD 
df$animals %>%
  filter(Outcome=="Succumbed") %>%
  select(Death) %>%
  unlist() %>%
  t.test()

# scores
# extract the challenge time and calculate time since challenge for all obs
baseline <- df$scores %>%
            group_by(NHP) %>%
            top_n(n=-1, wt=Datetime) %>%
            ungroup() %>%
            rename(Challenge=Datetime) %>%
            select(NHP, Challenge)
scores <- df$scores %>%
          # add in challenge dose
          left_join(df$animals, by="NHP") %>%
          # calculate time (days) since exact challenge time
          left_join(baseline, by="NHP") %>%
          mutate(DPI=as.integer(Datetime-Challenge),
                 DPI=DPI/(60*60*24)) %>%
          arrange(NHP)
scores %>%
  mutate(Score=if_else(Score > 9, 9, Score)) %>% # cap score at >=9 (euthanasia)
  ggplot(aes(DPI, Score)) +
  geom_hline(yintercept=9, col="lightgrey", linetype=3) +
  geom_line(aes(col=Group, group=NHP), lineend="round") +
  scale_color_manual(NULL, values=cols.group) +
  scale_x_continuous("Days postinfection", breaks=breaks.dates) +
  scale_y_continuous("Clinical score", breaks=c("0"=0, "3"=3, "6"=6, "≥9"=9)) +
  labs(title="Clinical scores") +
  guides(color=guide_legend(override.aes=list(linewidth=1))) +
  theme(legend.position=c(0.8, 0.6))
ggsave("analysis/score.png", units="in", width=3.75, height=2)

# who scored at least once? only succumbed and symptomatic
scores %>%
  filter(Score > 0) %>%
  select(NHP, Group) %>%
  distinct() %>%
  arrange(Group)

# when does scoring start?
scorestart <- scores %>%
              filter(Score > 0) %>%
              group_by(Group, NHP) %>%
              top_n(n=-1, wt=DPI) %>%
              ungroup()
# range and median by outcome
scorestart %>%
  group_by(Group) %>%
  summarise(Min=min(DPI),
            Max=max(DPI),
            Median=median(DPI))
# mean and 95% CI
scorestart %>%
  filter(Group=="Succumbed") %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
scorestart %>%
  filter(Group=="Symptomatic") %>%
  select(DPI) %>%
  unlist() %>%
  t.test()
# any difference between succumbed and symptomatic? No, mostly from high var.
scorestart %>%
  rstatix::wilcox_test(DPI ~ Outcome)
# clean up
rm(scores, scorestart)

# decreased appetite
appetite <- df$observations %>%
            select(NHP, DPI, `Biscuits eaten`) %>%
            left_join(df$animals, by="NHP") %>%
             # remove D0 since we don't count biscuits that day
            filter(DPI > 0)
# for baseline, take D1-D4 average since these NHPs had crazy variation
baseline <- appetite %>%
            filter(DPI %in% 1:4) %>%
            group_by(NHP) %>%
            summarise(Baseline=mean(`Biscuits eaten`),
                      .groups="drop")
# calculate appetite relative to baseline
appetite <- appetite %>%
            left_join(baseline, by="NHP") %>%
            mutate(Appetite=`Biscuits eaten`/Baseline) 
# max appetite capped at 1 for visualization purposes
appetite$Appetite[appetite$Appetite > 1] <- 1
# plot it -- need to shift DPI by -1 since we're counting biscuits from the 
# previous day
appetite %>%
  mutate(DPI=DPI-1) %>%
  ggplot(aes(DPI, NHP)) +
  geom_tile(aes(fill=Appetite), col="black") +
  scale_fill_gradient(NULL, low="black", high="white", na.value="black",
                      breaks=c(0, 0.5, 1), labels=c("0%", "50%", "100%")) +
  scale_x_continuous("Days postinfection", breaks=breaks.dates, 
                     limits=c(NA, 35), expand=c(0,0)) +
  scale_y_discrete("", expand=c(0,0)) +
  facet_wrap(~Group, ncol=1, scales="free_y") +
  labs(title="Appetite") +
  theme(legend.position=c(0.75, 0.9),
        legend.direction="horizontal",
        axis.text.y=element_blank(),
        axis.ticks.y=element_line(color="transparent"))
ggsave("analysis/appetite.png", units="in", width=3.75, height=4)

# clean up
rm(appetite, baseline)

## update LD50 -----------------------------------------------------------------
# using published data on this model https://doi.org/10.1016/j.xcrm.2025.102263

# load LD50 and make dose linear by taking log10
df$LD50 <- df$LD50 %>%
           mutate(Actual=log10(Actual),
                  Target=log10(Target))

# using target dose
drc::drm(Succumbed/Total ~ Target, weights=Total, data=df$LD50, 
         fct=drc::LL.2(names=c("Slope", "LD50")), 
         type="binomial", logDose=10) %>%
  drc::ED(c(10, 50, 90), display=FALSE) %>%
  round(digits=1)

# using back-titer dose
drc::drm(Succumbed/Total ~ Actual, weights=Total, data=df$LD50, 
         fct=drc::LL.2(names=c("Slope", "LD50")), 
         type="binomial", logDose=10) %>%
  drc::ED(c(10, 50, 90), display=FALSE) %>%
  round(digits=1)

## viremia ---------------------------------------------------------------------
# add challenge dose
viremia <- df$virology.blood %>%
           left_join(df$animals, by="NHP") %>%
           select(NHP, DPI, Outcome, `GEq/mL`, `PFU/mL`, Group) %>%
           arrange(Group)

# when was virus first detectable?
viremia %>%
  reshape2::melt(id.vars=c("NHP", "DPI", "Group", "Outcome")) %>%
  filter(value > 0) %>%
  group_by(NHP, variable) %>%
  top_n(n=-1, wt=DPI) %>%
  reshape2::dcast(NHP + Group + Outcome ~ variable, value.var="DPI") %>%
  arrange(Group, Outcome)

# qRT-PCR
viremia %>%
  mutate(`GEq/mL`=if_else(`GEq/mL`==0, 99, `GEq/mL`)) %>% # add 1e2 LOD to 0s
  ggplot(aes(DPI, `GEq/mL`)) +
  geom_line(aes(col=Group, group=NHP), lineend="round") +
  scale_color_manual(NULL, values=cols.group) +
  scale_y_continuous("N GEq/mL", limits=c(99, 1e13), 
                     breaks=c(99, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10") +
  scale_x_continuous("Days postinfection", breaks=breaks.dates) +
  labs(title="LASV RNA") +
  theme(legend.position="none")
ggsave("analysis/viremia-qpcr-longitudinal.png", 
       units="in", width=3.75, height=2)

# plaque assay
viremia %>%
  mutate(`PFU/mL`=if_else(`PFU/mL`==0, 24, `PFU/mL`)) %>% # add 24 LOD to 0s
  ggplot(aes(DPI, `PFU/mL`)) +
  geom_line(aes(group=NHP, col=Group), lineend="round") +
  scale_color_manual(NULL, values=cols.group) +
  scale_y_continuous("LASV PFU/mL", limits=c(24, 1e7), breaks=c(24, 1e3, 1e6),
                     labels=c("<LoD", "1e+03", "1e+06"), trans="log10") +
  scale_x_continuous("Days postinfection", breaks=breaks.dates) +
  labs(title="LASV titer") +
  theme(legend.position="none")
ggsave("analysis/viremia-plaqueassay-longitudinal.png",
       units="in", width=3.75, height=2)

# is there a difference in viral load at onset between symptomatic and succumbed
viremia %>%
  reshape2::melt(id.vars=c("NHP", "Outcome", "Group", "DPI"),
                 variable.name="Assay") %>%
  filter(value > 0,
         Group != "Asymptomatic") %>%
  group_by(NHP, Assay) %>%
  slice_min(n=1, order_by=DPI, with_ties=FALSE) %>%
  ungroup() %>%
  group_by(Assay) %>%
  rstatix::wilcox_test(value ~ Group)

# is there a difference in peak viral load?
# using all groups here since 0 is a valid value
viro.peaks <- viremia %>%
              select(-DPI) %>%
              reshape2::melt(id.vars=c("NHP", "Outcome", "Group"),
                             variable.name="Assay") %>%
              group_by(NHP, Assay) %>%
              slice_max(value, n=1, with_ties=FALSE) %>%
              ungroup() 
# check: what is min/max each group? Thresholding by outcome.
viro.peaks %>%
  group_by(Assay, Group) %>%
  summarise(Max=max(value),
            Min=min(value),
            .groups="drop")
# kw test first -- yes to both assays
viro.peaks %>%
  group_by(Assay) %>%
  rstatix::kruskal_test(value ~ Group)
# dunn post-hoc
pval <- viro.peaks %>%
        group_by(Assay) %>%
        rstatix::dunn_test(value ~ Group, p.adjust.method="none") %>%
        rstatix::add_xy_position(scales="free_y")
# plot these
viro.peaks %>%
  filter(Assay=="GEq/mL") %>%
  mutate(value=if_else(value==0, 99, value)) %>% # add 1e2 LOD to 0s
  ggplot(aes(Group, value)) +
  geom_hline(yintercept=5e10, linetype=3, col="grey", lineend="round") +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(fill=Group), height=0, pch=21) +
  scale_fill_manual(NULL, values=cols.group) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="GEq/mL"),
                             label="p.adj.signif", step.increase=0.12,
                             bracket.nudge.y=1) +
  scale_y_continuous("N GEq/mL", limits=c(NA, 1e16), 
                     breaks=c(99, 1e4, 1e8, 1e12),
                     labels=c("<LoD", "1e+04", "1e+08", "1e+12"),
                     trans="log10", expand=expansion(c(0.05, 0.1))) +
  labs(x=NULL, title="Max RNA") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
ggsave("analysis/viremia-qpcr-peak.png",
       units="in", width=3.75/2, height=2)
viro.peaks %>%
  filter(Assay=="PFU/mL") %>%
  mutate(value=if_else(value==0, 24, value)) %>% # add 24 LOD to 0s
  ggplot(aes(Group, value)) +
  geom_hline(yintercept=1e3, linetype=3, col="grey", lineend="round") +
  geom_boxplot(fill="grey", outliers=FALSE, alpha=0.5) +
  geom_jitter(aes(fill=Group), height=0, pch=21) +
  scale_fill_manual(NULL, values=cols.group) +
  ggpubr::stat_pvalue_manual(data=filter(pval, Assay=="PFU/mL"), 
                             label="p.adj.signif", step.increase=0.12,
                             bracket.nudge.y=1) +
  scale_y_continuous("LASV PFU/mL", limits=c(NA, 1e9), 
                     breaks=c(24, 1e3, 1e6),
                     labels=c("<LoD", "1e+03", "1e+08"),
                     trans="log10", expand=expansion(c(0.05, 0.1))) +
  labs(x=NULL, title="Max titer") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
ggsave("analysis/viremia-plaqueassay-peak.png",
       units="in", width=3.75/2, height=2)

# clean up
rm(viremia, viro.peaks, pval)

# tissue virology
df$virology.tissues %>%
  mutate(`GEq/g`=if_else(`GEq/g`==0, 1e4, `GEq/g`)) %>% # add LOD
  filter(Tissue != "Aqueous humor") %>%
  left_join(df$animals, by="NHP") %>%
  ggplot(aes(Tissue, `GEq/g`)) +
  geom_boxplot(aes(fill=Group, col=Group), alpha=0.5, na.rm=TRUE,
               outlier.size=0.5) +
  scale_fill_manual(NULL, values=cols.group) +
  scale_color_manual(NULL, values=cols.group) +
  scale_y_continuous("N GEq/mL", breaks=c(1e4, 1e8, 1e12, 1e16),
                     labels=c("<LoD", "1e+08", "1e+12", "1e+16"),
                     trans="log10", expand=expansion(c(0.1, 0.2))) +
  labs(x=NULL, title="LASV RNA in tissues") +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="bottom")
ggsave("analysis/viremia-qpcr-tissues.png",
       units="in", width=7.5, height=4)

## clinical pathology ----------------------------------------------------------
# format clinical chemistry
chem <- df$chemistry %>%
        left_join(select(df$animals, NHP, Outcome, Group, Death), 
                  by="NHP") %>%
        reshape2::melt(id.vars=c("NHP", "DPI", "Outcome", "Group", "Death"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        mutate(NHP=factor(NHP, levels=rev(levels(df$animals$NHP))),
               Type="chem")
# format hematology
hema <- df$hematology %>%
        select(-starts_with("Mean"), -starts_with("Hem")) %>%
        left_join(select(df$animals, NHP, Outcome, Group, Death),
                  by="NHP") %>%
        reshape2::melt(id.vars=c("NHP", "DPI", "Outcome", "Group", "Death"),
                       variable.name="Analyte",
                       value.name="Concentration") %>%
        mutate(Type="hema")
# combine and format
clinpath <- rbind(chem, hema)
chem <- levels(chem$Analyte)
hema <- levels(hema$Analyte)
# date a "datetime" variable for terminals
clinpath$Datetime <- clinpath$DPI
clinpath$Datetime[(clinpath$DPI != 35) & (clinpath$DPI==clinpath$Death)] <- "Terminal"
clinpath$Datetime <- factor(clinpath$Datetime, levels=unique(clinpath$Datetime))

# when do we see changes per analyte in succumbing animals?
# not all animals have distinct D14/terminal, so using KW instead of Friedman
# also, we can use all baseline instead of just n=4
analytes <- clinpath %>%
            filter(Group=="Succumbed" | Datetime=="0") %>%
            mutate(Datetime=droplevels(Datetime)) %>%
            group_by(Analyte) %>%
            rstatix::kruskal_test(Concentration ~ Datetime) %>%
            filter(p < 0.05) 
# run Dunn post-hoc and filter for comparisons relative to D0/baseline
clinpath %>%
  filter(Group=="Succumbed" | Datetime=="0",
         Analyte %in% analytes$Analyte) %>%
  mutate(Datetime=droplevels(Datetime)) %>%
  group_by(Analyte) %>%
  rstatix::dunn_test(Concentration ~ Datetime) %>%
  select(Analyte, group1, group2, p.adj, p.adj.signif) %>%
  rstatix::p_format(new.col=TRUE, digits=1) %>%
  filter(group1==0,
         p.adj < 0.05)

# what about symptomatic animals?
# we could use friedman, but need to compare to succumbed AGMs -- 
# using KW to be consistent
analytes <- clinpath %>%
            filter(Group=="Symptomatic" | Datetime=="0") %>%
            mutate(Datetime=droplevels(Datetime)) %>%
            group_by(Analyte) %>%
            rstatix::kruskal_test(Concentration ~ Datetime) %>%
            filter(p < 0.05) 
# run Dunn post-hoc and filter for comparisons relative to D0/baseline
# ALB is nearly significant at D14 but we have only 6x samples (3:3)
clinpath %>%
  filter(Group=="Symptomatic" | Datetime=="0",
         Analyte %in% analytes$Analyte) %>%
  mutate(Datetime=droplevels(Datetime)) %>%
  group_by(Analyte) %>%
  rstatix::dunn_test(Concentration ~ Datetime) %>%
  select(Analyte, group1, group2, p.adj, p.adj.signif) %>%
  rstatix::p_format(new.col=TRUE, digits=1) %>%
  filter(group1==0,
         p.adj < 0.05)

# repeat for asymptomatic
analytes <- clinpath %>%
            filter(Group=="Asymptomatic" | Datetime=="0") %>%
            mutate(Datetime=droplevels(Datetime)) %>%
            group_by(Analyte) %>%
            rstatix::kruskal_test(Concentration ~ Datetime) %>%
            filter(p < 0.05) 
# run Dunn post-hoc and filter for comparisons relative to D0/baseline
# ALB is nearly significant at D14 but we have only 6x samples (3:3)
clinpath %>%
  filter(Group=="Asymptomatic" | Datetime=="0",
         Analyte %in% analytes$Analyte) %>%
  mutate(Datetime=droplevels(Datetime)) %>%
  group_by(Analyte) %>%
  rstatix::dunn_test(Concentration ~ Datetime) %>%
  select(Analyte, group1, group2, p.adj, p.adj.signif) %>%
  rstatix::p_format(new.col=TRUE, digits=1) %>%
  filter(group1==0,
         p.adj < 0.05)
rm(analytes)

# plot all spaghetti
spaghetti <- c(chem, hema) %>%
             lapply(function(i) {
               # plot it
               clinpath %>%
                 filter(Analyte==i,
                        !is.na(Concentration)) %>%
                 ggplot(aes(DPI, Concentration)) +
                 geom_line(aes(group=NHP, col=Group), lineend="round") +
                 scale_color_manual(NULL, values=cols.group) +
                 scale_x_continuous("Days postinfection", breaks=breaks.dates) +
                 labs(y=df$units$Units[df$units$Analyte==i], title=i) +
                 theme(legend.position="none")
            })
names(spaghetti) <- c(chem, hema)
# plot clinical chemistry; pull out legend and paste on the bottom
legs <- spaghetti$GLU +
        guides(color=guide_legend(override.aes=list(linewidth=1))) +
        theme(legend.position="bottom")
legs <- cowplot::get_legend(legs)
legs <- cowplot::ggdraw(legs)
x <- cowplot::plot_grid(plotlist=spaghetti[chem], 
                        labels="AUTO", ncol=3)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(17, 1))
ggsave("analysis/clinpath-spaghetti-chemistry.png",
       units="in", width=7.5, height=7)
# repeat for hematology; paste legend in the bottom right corner
legs <- spaghetti$GLU +
        guides(color=guide_legend(title=NULL, override.aes=list(linewidth=1))) +
        theme(legend.position="right")
legs <- cowplot::get_legend(legs)
spaghetti$legs <- cowplot::ggdraw(legs)
cowplot::plot_grid(plotlist=spaghetti[c(hema, "legs")], 
                   labels=LETTERS[1:length(hema)], ncol=3)
ggsave("analysis/clinpath-spaghetti-hematology.png",
       units="in", width=7.5, height=6)
# plot selected spaghetti
cowplot::plot_grid(plotlist=spaghetti[c("ALB", "Platelets", "Neutrophils", 
                                        "CRP", "AST", "legs")],
                   nrow=2, labels=LETTERS[1:5])
ggsave("analysis/clinpath-spaghetti-selected.png",
       units="in", width=7.5, height=4)
rm(legs, x, spaghetti)

# calculate and compare baseline (pooled), max, and min
# baseline
comparison <- clinpath %>%
              filter(DPI==0) %>%
              mutate(Comparison="D0")
# max
comparison <- clinpath %>%
              group_by(NHP, Analyte) %>%
              slice_max(order_by=Concentration, with_ties=FALSE) %>%
              mutate(Comparison="Max") %>% 
              rbind(comparison)
# min
comparison <- clinpath %>%
              group_by(NHP, Analyte) %>%
              slice_min(order_by=Concentration, with_ties=FALSE) %>%
              mutate(Comparison="Min") %>% 
              rbind(comparison) %>%
              mutate(Comparison=factor(Comparison, 
                                       levels=c("D0", "Min", "Max")),
                     # update group D0
                     Group=factor(Group, levels=c("D0", names(cols.group))),
                     Group=if_else(Comparison=="D0", "D0", Group),
                     Group=factor(Group, levels=c("D0", names(cols.group)))) %>%
              ungroup()

# insufficient samples (<8) for symptomatic vs. succumbed comparisons; 
# use thresholding
comparison %>%
  filter(Group!="Asymptomatic") %>%
  group_by(Analyte, Group, Comparison) %>%
  summarise(Max=max(Concentration),
            Min=min(Concentration),
            .groups="drop") %>%
  reshape2::melt(id.vars=c("Group", "Analyte", "Comparison"),
                 value.name="Concentration") %>%
  mutate(variable=interaction(Group, variable),
         Comparison=if_else(Comparison=="D0", 
                            str_remove(variable, "D0\\."),
                            Comparison)) %>%
  reshape2::dcast(Analyte + Comparison ~ variable, value.var="Concentration") %>%
  filter((Comparison=="Max" &
            (Succumbed.Min > Symptomatic.Max | 
             Symptomatic.Min > Succumbed.Max)) |
         (Comparison=="Min" & 
            (Succumbed.Max < Symptomatic.Min |
             Symptomatic.Max < Succumbed.Min)))

# comparison #1: max/min each group vs. pooled baseline (MWU/WX)
# assumptions check: only keep max/min groups with >1 value (not all <LOD)
comps.bl <- comparison %>%
            filter(Comparison != "D0") %>%
            select(Analyte, Comparison, Group, Concentration) %>%
            distinct() %>%
            group_by(Analyte, Comparison, Group) %>%
            summarise(Values=n(),
                      .groups="drop") %>%
            filter(Values > 1) %>%
            mutate(groups=interaction(Analyte, Comparison, Group)) %>%
            select(groups) %>%
            unlist() %>%
            as.character()
comps.bl <- comparison %>%
            mutate(groups=interaction(Analyte, Comparison, Group)) %>%
            filter(Comparison=="D0" | groups %in% comps.bl) %>%
            mutate(groups=interaction(Comparison, Group)) %>% 
            group_by(Analyte) %>%
            rstatix::wilcox_test(Concentration ~ groups, ref.group="D0.D0") %>%
            filter(p.adj < 0.05) %>%
            # add x positions
            mutate(xmin=1, # D0
                   # extract comparison and group
                   Comparison=factor(str_extract(group2, "^[^\\.]+"),
                                     levels=levels(comparison$Comparison)),
                   Group=factor(str_extract(group2, "[^\\.]+$"),
                                levels=names(cols.group)),
                   # format xmax position
                   xmax=factor(Group, labels=c(-0.266667, 0, 0.266667)),
                   xmax=as.numeric(as.character(xmax)),
                   xmax=as.numeric(Comparison)+xmax) %>%
            # arrange by the rightmost x position
            arrange(xmax)

# comparison #2: inter-group comparison at max/min (KW -> dunn)
comps.grps <- comparison %>%
              filter(Group != "D0") %>%
              group_by(Analyte, Comparison) %>%
              rstatix::kruskal_test(Concentration ~ Group) %>%
              # keep only valid comparisons
              filter(p < 0.05) %>%
              select(Analyte, Comparison) %>%
              # add back data to run dunn
              left_join(comparison, by=c("Analyte", "Comparison")) %>%
              group_by(Analyte, Comparison) %>%
              rstatix::dunn_test(Concentration ~ Group) %>%
              filter(p.adj < 0.05) %>%
              # add x positions
              rstatix::add_xy_position(x="Comparison", step.increase=0)

# plot everything
comparison <- c(chem, hema) %>%
              lapply(function(i) {
                # subset comparisons and stats
                analyte <- filter(comparison, Analyte==i) # data
                bl <- filter(comps.bl, Analyte==i) # comparisons to baseline
                grp <- filter(comps.grps, Analyte==i) # comparisons within peaks
                
                # build the basic plot
                p <- analyte %>%
                     ggplot(aes(Comparison, Concentration)) +
                     geom_boxplot(aes(group=interaction(Group, Comparison)), 
                                  fill="lightgrey", alpha=0.5, col="black", 
                                  outliers=FALSE) +
                     geom_point(aes(group=interaction(Group, Comparison), 
                                    col=Group), size=0.5, 
                                position=position_jitterdodge(jitter.width=0.6))
                
                # add comparisons to baseline, if any
                if(dim(bl)[1] > 0) {
                  # define increment and add y coordinates
                  increment <- 0.2*max(analyte$Concentration)
                  bl$y.position <- max(analyte$Concentration) + 
                                   increment*1:dim(bl)[1]  
                  # add to plot
                  p <- p +
                       ggpubr::stat_pvalue_manual(data=bl, color="Group", 
                                                  tip.length=0)
                }
                
                # add inter-group comparisons, if any
                if(dim(grp)[1] > 0) {
                  p <- p +
                       ggpubr::stat_pvalue_manual(data=grp, step.increase=0.2, 
                                                  step.group.by="Comparison", 
                                                  bracket.nudge.y=0.5*increment)
                }
                
                # finish building plot and return
                p +
                  scale_color_manual(values=cols.group) +
                  scale_y_continuous(df$units$Units[df$units$Analyte==i],
                                     expand=expansion(mult=c(0.05, 0.1))) +
                  labs(x=NULL,
                       title=i) +
                  theme(legend.position="none")
              })
names(comparison) <- c(chem, hema)
# plot clinical chemistry; pull out legend and paste on the bottom
legs <- comparison$Erythrocytes +
        guides(color=guide_legend(title=NULL, override.aes=list(size=2))) +
        theme(legend.position="bottom")
legs <- cowplot::get_legend(legs)
legs <- cowplot::ggdraw(legs)
x <- cowplot::plot_grid(plotlist=comparison[chem], 
                        labels="AUTO", ncol=4)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(17, 1))
ggsave("analysis/clinpath-peak-chemistry.png",
       units="in", width=7.5, height=7)
# repeat for hematology; paste legend in the bottom right corner
x <- cowplot::plot_grid(plotlist=comparison[hema], 
                        labels="AUTO", ncol=4)
cowplot::plot_grid(x, legs, ncol=1, rel_heights=c(12, 1))
#cowplot::plot_grid(plotlist=comparison[c(hema, "legs")], 
#                   labels=LETTERS[1:length(hema)], ncol=4)
ggsave("analysis/clinpath-peak-hematology.png",
       units="in", width=7.5, height=5)
# plot selected spaghetti
cowplot::plot_grid(plotlist=comparison[c("ALB", "Platelets", "Neutrophils", 
                                         "CRP", "AST", "legs")],
                   nrow=2, labels=LETTERS[6:10])
ggsave("analysis/clinpath-peak-selected.png",
       units="in", width=7.5, height=4)

# clean up
rm(legs, x, comparison, comps.bl, comps.grps, chem, hema, clinpath)

## ELISAs ----------------------------------------------------------------------
# over time -- very messy
df$ELISAs %>%
  reshape2::melt(measure.vars=c("IgM", "IgG"),
                 value.name="Titer") %>%
  left_join(df$animals, by="NHP") %>%
  ggplot(aes(DPI, Titer+100)) +
  geom_line(aes(group=NHP, col=Group)) +
  geom_point(aes(fill=Group), pch=21) +
  scale_fill_manual(values=cols.group) +
  scale_color_manual(values=cols.group) +
  scale_y_continuous(limits=c(NA, 30000),
                     breaks=c(100, 200, 800, 3200, 12800),
                     labels=c("<LOD", 200, 800, 3200, 12800),
                     trans="log2") +
  facet_grid(variable ~ Antigen)

# is there any difference between groups at peak titers?
elisas <- df$ELISAs %>%
          group_by(NHP, Antigen) %>%
          top_n(n=1, wt=DPI) %>%
          reshape2::melt(measure.vars=c("IgM", "IgG"),
                         variable.name="Immunoglobulin",
                         value.name="Titer") %>%
          left_join(df$animals, by="NHP")

# get p-value significance by testing by KW and then dunn post-hoc
elisas %>%
  mutate(Titer=log2(Titer+1)) %>%
  group_by(Immunoglobulin, Antigen) %>%
  rstatix::kruskal_test(Titer ~ Group) %>%
  filter(p < 0.05)
# add pseudocount to samples <LOD
elisas$Titer[elisas$Titer==0] <- 50
# plot it and add statistical comparisons
elisas %>%
  ggplot(aes(Antigen, Titer)) +
  geom_boxplot(aes(group=interaction(Antigen, Group)), 
               fill="lightgrey", col="black",
               alpha=0.5, outliers=FALSE) +
  geom_point(aes(group=interaction(Antigen, Group), size=Group, col=Group), 
             position=position_jitterdodge(), size=0.5) +
  scale_color_manual(NULL, values=cols.group) +
  ggpubr::stat_pwc(aes(group=Group), method="dunn_test", hide.ns=TRUE, 
                   p.adjust.by="panel", label="p.adj.signif") +
  scale_y_continuous("Endpoint titer", breaks=c(50, 200, 800, 3200, 12800),
                     labels=c("<LOD", 200, 800, 3200, 12800),
                     trans="log2", expand=expansion(c(0.05, 0.2))) +
  facet_wrap(~Immunoglobulin) +
  labs(x="LASV antigen", title="ELISAs") +
  theme(legend.position="none")
ggsave("analysis/elisas-endpoint.png",
       units="in", width=3, height=2.5)

# clean up
rm(psig, elisas)

## PRNTs -----------------------------------------------------------------------
# get the average plaques per sample
prnt <- df$PRNTs %>%
        reshape2::melt(measure.vars=c("Plaques A", "Plaques B"),
                       value.name="Plaques") %>%
        group_by(NHP, DPI, Batch, `Dilution factor`) %>%
        summarise(Plaques=mean(Plaques),
                  .groups="drop")
# pull out VC and calculate the baseline for each batch
vc <- prnt %>%
      filter(NHP=="VC") %>%
      rename(VC=Plaques) %>%
      select(Batch, VC)
# calculate % reduction
prnt <- prnt %>%
        filter(NHP != "VC") %>%
        # add in VC counts and calculate % reduction
        left_join(vc, by="Batch") %>%
        # add NHP info
        left_join(df$animals, by="NHP") %>%
        # calculate % reduction
        mutate(Reduction=VC-Plaques,
               PercentReduction=100*Reduction/VC)
# update the DPI to factor (0, 14, "terminal", 35)
prnt$DPI[(prnt$Outcome=="Succumbed") & (prnt$DPI==prnt$Death)] <- "Terminal"
prnt$DPI <- factor(prnt$DPI, levels=c("0", "14", "Terminal", "35"))

# bound the % reduction to [0, 100]
prnt$PercentReduction[prnt$PercentReduction < 0] <- 0

# plot each NHP by group
prnt %>%
  arrange(desc(DPI)) %>%
  ggplot(aes(`Dilution factor`, PercentReduction)) +
  geom_hline(yintercept=50, linetype=3, col="lightgrey") +
  geom_line(aes(group=interaction(NHP, DPI), col=DPI)) +
  geom_point(aes(fill=DPI), col="black", size=2, pch=21) +
  scale_fill_manual(values=c("white", "darkgrey", 
                             "#e41a1c", "black")) +
  scale_color_manual(values=c("white", "darkgrey", 
                              "#e41a1c", "black")) +
  scale_x_continuous(trans="log2") +
  scale_y_continuous("Neutralization (%)", limits=c(0, 100)) +
  facet_wrap(~Group, nrow=1) +
  labs(title="PRNT curves") +
  guides(fill=guide_legend(override.aes=list(pch=21)),
         shape=guide_legend(override.aes=list(fill="black")))
ggsave("analysis/prnts-curves.png",
       units="in", width=7.5, height=3)

# plot neutralization for each at the lowest dilution and compare outcomes
prnt <- prnt %>%
        filter(`Dilution factor`==10,
               DPI %in% c("Terminal", "35")) %>%
        group_by(NHP) %>%
        slice_max(order_by=PercentReduction, n=1, with_ties=FALSE) %>%
        ungroup()
# all groups by KW to check for significance
prnt %>%
  rstatix::kruskal_test(PercentReduction ~ Group)
# dunn's post hoc
pval <- prnt %>%
        rstatix::dunn_test(PercentReduction ~ Group) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x="Group")
prnt %>%
  ggplot(aes(Group, PercentReduction)) +
  geom_hline(yintercept=50, linetype=3, col="lightgrey") +
  geom_boxplot(alpha=0.5, outliers=FALSE, col="black", fill="lightgrey") +
  geom_jitter(aes(fill=Group), height=0, pch=21) +
  scale_fill_manual(NULL, values=cols.group) +
  scale_y_continuous("Neutralization (%)", breaks=c(0, 25, 50, 75, 100),
                     expand=expansion(mult=c(0.1, 0.2))) +
  ggpubr::stat_pwc(aes(group=Group), method="dunn_test", hide.ns=TRUE, 
                   label="p.adj.signif") +
  labs(x=NULL, title="LASV neutralization") +
  guides(fill=guide_legend(override.aes=list(size=3))) +
  theme(axis.text.x=element_text(angle=15, hjust=0.6, vjust=0.7))
ggsave("analysis/prnts-max.png",
       units="in", width=4.5, height=2.5)

# clean up
rm(prnt, vc)

## temperature loggers ---------------------------------------------------------
logs <- df$loggers %>%
        # extract date from Datetime and update to DPI
        mutate(Datetime=as.Date(Datetime),
               DPI=difftime(Datetime, "2025-03-19", units="days"),
               DPI=round(as.numeric(DPI))) %>%
        # get the mean temperature and activity score per day 
        group_by(NHP, DPI) %>%
        summarise(Temperature=mean(Temperature, na.rm=TRUE),
                  Activity=mean(ACT, na.rm=TRUE),
                  .groups="drop") %>%
        # add metadata
        left_join(df$animals, by="NHP")

# calculate baseline distribution (mean and StDev of prechallenge means)
tempdist <- logs %>%
            filter(DPI < 0) %>%
            summarise(Mean=mean(Temperature),
                      StDev=sd(Temperature)) %>%
            unlist()
tempdist <- data.frame(DPI=0:35,
                       Max=tempdist["Mean"] + 2*tempdist["StDev"],
                       Min=tempdist["Mean"] - 2*tempdist["StDev"])

# plot mean daily temperature
logs %>%
  filter(DPI >= 0) %>%
  ggplot(aes(DPI)) +
  # draw ribbon showing 95% CI (mean +/- 2x stdev)
  geom_ribbon(data=tempdist, aes(ymin=Min, ymax=Max), 
              fill="lightgrey", alpha=0.5) +
  geom_line(aes(y=Temperature, group=interaction(NHP, Group), col=Group), 
            lineend="round") +
  scale_color_manual(values=cols.group) +
  scale_x_continuous("Days postinfection", breaks=breaks.dates) +
  scale_y_continuous("Temperature (°C)", limits=c(33, NA)) +
  labs(col=NULL,
       title="Mean daily temperature") +
  theme(legend.position=c(0.7, 0.25))
ggsave("analysis/temperature.png",
       units="in", width=3.75, height=2.5)

# clean up
rm(logs, tempdist)

## fin -------------------------------------------------------------------------
sessionInfo()
