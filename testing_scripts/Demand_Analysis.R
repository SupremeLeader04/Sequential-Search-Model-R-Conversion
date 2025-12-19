total_purchase = 0
total_searched = 0
total_searches = 0
source('simWeitz.R')
source('simWeitzorig.R')

for (i in 1:50) {
  print(i)
  data = simWeitz(1000, 5, c(-1.0464921 ,  1.46654103,  1.11547285, -1.26343864,  0.24762589), i)
  data_purchase = data[data[, 2] == 1, ncol(data)]
  total_purchase = total_purchase + sum(data_purchase) / length(data_purchase)
  data_searched = data[, ncol(data) - 3]
  total_searched = total_searched + sum(data_searched) / length(data_searched)
  data_searches = data[, ncol(data) - 1]
  total_searches = total_searches + sum(data_searches) / 1000
}

print('Average Proportion Demanding:')
print(1 - (total_purchase / 50))
print('Average Proportion Searching at Least Once:')
print(total_searched / 50)
print('Average Number of Searches per Consumer:')
print(total_searches / 50)