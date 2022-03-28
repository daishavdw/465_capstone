library(tidyverse)
library(ggvenn)

clinGenData = read_csv("clinGen.csv")
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

caddGenes = scan("caddGenes", what="")
numCadd = read_tsv("numCadd")
caddScores = read_tsv("caddScores")

sourceNum = bind_rows(numClinGen, numSift)
sourceNum = bind_rows(sourceNum, numCadd)


ggplot(data=sourceNum) +
  geom_col(aes(x=Source, y=Count)) +
  geom_text(aes(x=Source, y=Count, label = Count), vjust = -0.5) +
  theme_bw()

geneVenn = list(ClinGen = clinGenGenes, SIFT = siftGenes, CADD = caddGenes)

ggvenn(data=geneVenn, show_percentage = FALSE)

ggplot() +
  geom_col(data=fig3Clin, aes(x=Assertion, y=Num)) +
  labs(y = "Count") +
  scale_x_discrete(labels = function(Assertion) str_wrap(Assertion, width = 10)) +
  theme_bw()

ggplot(data=siftScores, aes(x=factor(Value), y=Count)) +
  geom_col(position = position_nudge(x = -0.5), width=1) +
  labs(x = "SIFT Score") +
  theme_bw()

ggplot(data=caddScores, aes(x=factor(Score), y=Count)) +
  geom_col(position = position_nudge(x = -0.5), width=1) +
  labs(x = "Phred Score") +
  theme_bw()
