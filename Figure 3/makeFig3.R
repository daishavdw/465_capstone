library(tidyverse)
library(ggvenn)

clinGenData = read_csv("clinGenPathogenicity.csv")
clinGen = filter(clinGenData, parseMethod != "ERROR" & pos != "ISSUE")

numClinGen = summarize(clinGen, Count = n()) %>%
  mutate(Source = "ClinGen") %>%
  relocate(Source, .before = Count)

clinGenGenes = unique(pull(clinGen, "HGNC Gene Symbol"))

fig3Clin = select(clinGen, Assertion) %>%
  group_by(Assertion) %>%
  summarize(Num = n())

siftGenes = scan("siftGenes", what = "")
numSift = read_tsv("numSift")
siftScores = read_tsv("siftScoreCounts")

sourceNum = bind_rows(numClinGen, numSift)



ggplot(data=sourceNum) +
  geom_col(aes(x=Source, y=Count)) +
  geom_text(aes(x=Source, y=Count, label = Count), vjust = -0.5) +
  theme_bw()

geneVenn = list(ClinGen = clinGenGenes, SIFT = siftGenes)

ggvenn(data=geneVenn)

ggplot() +
  geom_col(data=fig3Clin, aes(x=Assertion, y=Num)) +
  labs(y = "Count") +
  scale_x_discrete(labels = function(Assertion) str_wrap(Assertion, width = 10)) +
  theme_bw()

ggplot(data=siftScores, aes(x=factor(Value), y=Count)) +
  geom_col(position = position_nudge(x = -0.5), width=1) +
  labs(x = "SIFT Score") +
  theme_bw()
  


