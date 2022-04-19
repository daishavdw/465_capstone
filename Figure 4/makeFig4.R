library(tidyverse)
library(ggvenn)

clinGenData = read_csv("clinGen.csv")
clinGen = filter(clinGenData, parseMethod != "ERROR" & pos != "ISSUE")

numClinGen = summarize(clinGen, Count = n()) %>%
  mutate(Source = "ClinGen") %>%
  relocate(Source, .before = Count)

clinGenGenes = unique(pull(clinGen, "HGNC Gene Symbol"))
clinGenGenes <- clinGenGenes[!clinGenGenes %in% c("N/A", "-")]

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


ggplot(data=sourceNum, aes(x=Source, y=Count, label = format(Count, big.mark = ",", scientific = FALSE))) +
  geom_col() +
  geom_text(nudge_x = -0.05, vjust= -0.5) +
  theme_bw()

geneVenn = list(ClinGen = clinGenGenes, SIFT = siftGenes, CADD = caddGenes)

ggvenn(data=geneVenn, show_percentage = FALSE)

ggplot() +
  geom_col(data=fig3Clin, aes(x=Assertion, y=Num)) +
  labs(y = "Count", x = "ClinGen Assertion") +
  scale_x_discrete(labels = function(Assertion) str_wrap(Assertion, width = 10)) +
  theme_bw()

ggplot(data=siftScores, aes(x=factor(Value), y=Count)) +
  geom_col(position = position_nudge(x = -0.5), width=1) +
  labs(x = "SIFT Score") +
  theme_bw()

ggplot(data=caddScores, aes(x=factor(Score), y=log(Count))) +
  geom_col(position = position_nudge(x = -0.5), width=1) +
  labs(x = "CADD Phred-like C Score", y = "log(Count)") +
  theme_bw()
