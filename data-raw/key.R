d2r <- readRDS("../data-clean/d2r.rds")

key <- merge(merge(aggregate(grp.id ~ grp.state, data=d2r, unique),
aggregate(grp.region ~ grp.state, data=d2r, unique)),
aggregate(grp.division ~ grp.state, data=d2r, unique))

key$id <- key$grp.id; key$grp.id <- NULL
key$state <- key$grp.state; key$grp.state <- NULL
key$division <- key$grp.division; key$grp.division <- NULL
key$region <- key$grp.region; key$grp.region <- NULL

fips <- c("01","02","04","05","06","08","09","10","12","13","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","44","45","46","47","48","49","50","51","53","54","55","56","60","66","69","72","78","11")
state <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming", "American Samoa", "Guam", "Northern Mariana Islands", "Puerto Rico", "Virgin Islands", "District of Columbia")
df <- data.frame(state,fips)

key <- merge(key,df)
key <- key[order(key$id),]
key

saveRDS(key,file="../data-clean/key.rds")
