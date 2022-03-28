library(tidyverse)

siftScores = read_tsv("siftScores")

table = filter(siftScores, SIFT_score <= .100) %>%
  summarize("0.100" = n()) %>%
  print()

hundreths2 = filter(siftScores, SIFT_score > .100 & SIFT_score <= .200) %>%
  summarize("0.200" = n()) %>%
  print()

table = bind_cols(table, hundreths2)

hundreths3 = filter(siftScores, SIFT_score > .200 & SIFT_score <= .300) %>%
  summarize("0.300" = n()) %>%
  print()

table = bind_cols(table, hundreths3)

hundreths4 = filter(siftScores, SIFT_score > .300 & SIFT_score <= .400) %>%
  summarize("0.400" = n()) %>%
  print()

table = bind_cols(table, hundreths4)

hundreths5 = filter(siftScores, SIFT_score > .400 & SIFT_score <= .500) %>%
  summarize("0.500" = n()) %>%
  print()

table = bind_cols(table, hundreths5)

hundreths6 = filter(siftScores, SIFT_score > .500 & SIFT_score <= .600) %>%
  summarize("0.600" = n()) %>%
  print()

table = bind_cols(table, hundreths6)

hundreths7 = filter(siftScores, SIFT_score > .600 & SIFT_score <= .700) %>%
  summarize("0.700" = n()) %>%
  print()

table = bind_cols(table, hundreths7)

hundreths8 = filter(siftScores, SIFT_score > .700 & SIFT_score <= .800) %>%
  summarize("0.800" = n()) %>%
  print()

table = bind_cols(table, hundreths8)

hundreths9 = filter(siftScores, SIFT_score > .800 & SIFT_score <= .900) %>%
  summarize("0.900" = n()) %>%
  print()

table = bind_cols(table, hundreths9)

hundreths10 = filter(siftScores, SIFT_score > .900 & SIFT_score <= 1.000) %>%
  summarize("1.000" = n()) %>%
  print()

table = bind_cols(table, hundreths10)

print(table)

table = pivot_longer(table, cols=1:10, names_to = "Value", values_to = "Count")

write_tsv(table, "siftScoreCounts")
