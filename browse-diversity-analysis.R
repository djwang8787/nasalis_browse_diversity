pacman::p_load(dplyr, ggplot2, randomForest, lme4, lubridate, MuMIn, vegan, reshape2, reshape, performance,
               see)

#### 1.1 data prep ####
data = read.csv('browse-data.csv')

data = data %>%
  mutate_at(c('Browse.Species', "ID", "Weather", "Leaves", "Avoid"), factor) %>%
  mutate(Date = as.POSIXct(data$Date, format = "%d/%m/%y"),
         Consumed = as.numeric(data$Consumed))

# H = shannons, D = Simpsons
data = data %>%
  filter(!ID %in% c("")) %>%
  group_by(Session, ID) %>%
  mutate(H = vegan::diversity(Consumed, index = "shannon"),
         D = vegan::diversity(Consumed, index = "simpson"))

data = data %>%
  filter(!ID %in% c("")) %>%
  group_by(ID, Session) %>%
  mutate(total = sum(Consumed)) 

data = data %>%
  group_by(ID, Session) %>%
  rowwise() %>%
  mutate(prop = (Consumed / total))

# fnb = food niche breadth
data = data %>%
  group_by(ID, Session) %>%
  mutate(
    sum.prop = sum(prop^2),
         Pi = (sum.prop^-1)-1, 
         N = n(),
         fnb = Pi/(N-1))

#### Monkey data ####
pm.overview = data %>%
  filter(Browse.Species == "Ketapung") %>%
  group_by(ID) %>%
  summarise(Shannons = mean(H),
            Simpsons = mean(D),
            FNB = mean(fnb),
            Ketapung = sum(Consumed))
pm2 = data %>%
  group_by(ID) %>%
  summarise(total = sum(Consumed))

pm.overview$Total = pm2$total
pm.overview$Age = c(3, 13, 3, 14, 23, 12, 12)
pm.overview$Sex = c("F", "M", "F", "F", "F", "F", "F")


#### Browse distribution ####
browse.d = data %>%
  filter(Sequence == 1) %>%
  group_by(Browse.Species) %>%
  summarise(availability = n())

browse.d = browse.d %>%
  filter(!Browse.Species == "NV") %>%
  rowwise() %>%
  mutate(prop = (availability/39)*100)
  
browse.d %>%
  ggplot(., aes(x = reorder(Browse.Species, prop), y = prop)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  xlab("Browse species") +
  ylab("Browse availability (%)")+
  coord_flip() 

#### consumption distribution ####
consumed.d = data %>%
  group_by(Browse.Species) %>%
  summarise(consumed = sum(Consumed))

consumed.d = consumed.d %>%
  mutate(total = sum(consumed)) %>%
  group_by(Browse.Species) %>%
  mutate(consume.prop = consumed / total * 100)

consumed.d %>%
  ggplot(., aes(x = reorder(Browse.Species, consume.prop), y = consume.prop)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  xlab("Browse species") +
  ylab("Browse consumed (%)")+
  coord_flip() 


#### 1.3 Random Forest models ####
rf.data = data %>%
  dcast(ID + H + Session + Leaves + Weather + fnb + D + total ~ Browse.Species, 
       value.var = "Consumed",
       fun.aggregate = function(x) if(length(x) == 0) NA_real_ else sum(x, na.rm = TRUE))


#### Shannon's Index ####
h.data = rf.data %>%
  select(-NV, -D, -fnb, -total)

h.data = rfImpute(H~., h.data)

H.model = h.data %>%
  randomForest(H ~., 
               data = .,
               mtry = sqrt(nrow(.)),
               method = "regression", 
               ntree = 500,
               importance = TRUE)

H.model
# Residuals are 0.18; H indexs are typically predicted wrong by 0.19

H.ImpData <- as.data.frame(randomForest::importance(H.model))
H.ImpData$Var.Names <- row.names(H.ImpData)

ggplot(H.ImpData, aes(x=reorder(Var.Names, `%IncMSE`), y=`%IncMSE`)) +
  geom_segment( aes(x=reorder(Var.Names, `%IncMSE`), xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

#### Simpson's Index ####

d.data = rf.data %>%
  select(-NV, -H, -fnb, -total)

d.data = rfImpute(D~., d.data)

D.model = d.data %>%
  randomForest(D ~., 
               data = .,
               mtry = sqrt(nrow(.)),
               method = "regression", 
               ntree = 500,
               importance = TRUE)


D.model
# Residuals are 0.05; H index are typically predicted wrong by 0.04

D.ImpData <- as.data.frame(randomForest::importance(D.model))
D.ImpData$Var.Names <- row.names(D.ImpData)

ggplot(D.ImpData, aes(x=reorder(Var.Names, `%IncMSE`), y=`%IncMSE`)) +
  geom_segment( aes(x=reorder(Var.Names, `%IncMSE`), xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

#### FNB index ####
fnb.data = rf.data %>%
  select(-NV, -H, -D, -total)

fnb.data = rfImpute(fnb~., fnb.data)

fnb.model = fnb.data %>%
  randomForest(fnb ~., 
               data = .,
               mtry = sqrt(nrow(.)),
               method = "regression", 
               ntree = 500,
               importance = TRUE)


fnb.model
# Residuals are 0.01; FNB index are typically predicted wrong by 0.009

FNB.ImpData <- as.data.frame(randomForest::importance(fnb.model))
FNB.ImpData$Var.Names <- row.names(FNB.ImpData)

ggplot(FNB.ImpData, aes(x=reorder(Var.Names, `%IncMSE`), y=`%IncMSE`)) +
  geom_segment( aes(x=reorder(Var.Names, `%IncMSE`), xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

#### Total consumption ####
tc.data = rf.data %>%
  select(-NV, -H, -D, -fnb)

tc.data = rfImpute(total~., tc.data)

tc.model = tc.data %>%
  randomForest(total ~., 
               data = .,
               mtry = sqrt(nrow(.)),
               method = "regression", 
               ntree = 500,
               importance = TRUE)


mf.tcdata = missForest::missForest(tc.data)
# Residuals are 117

tc.ImpData <- as.data.frame(randomForest::importance(tc.model))
tc.ImpData$Var.Names <- row.names(tc.ImpData)

ggplot(tc.ImpData, aes(x=reorder(Var.Names, `%IncMSE`), y=`%IncMSE`)) +
  geom_segment( aes(x=reorder(Var.Names, `%IncMSE`), xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

#### KTP consumption ####
ktp.data = rf.data %>%
  select(-NV, -H, -D, -fnb, -total)

ktp.data = rfImpute(Ketapung~., ktp.data)

ktp.model = ktp.data %>%
  randomForest(Ketapung ~., 
               data = .,
               mtry = sqrt(nrow(.)),
               method = "regression", 
               ntree = 500,
               importance = TRUE)


ktp.model
# Residuals are 57

ktp.ImpData <- as.data.frame(randomForest::importance(ktp.model))
ktp.ImpData$Var.Names <- row.names(ktp.ImpData)

ggplot(ktp.ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# a more beautiful plot!
ktp.ImpData %>%
ggplot(., aes(x=reorder(Var.Names, `%IncMSE`), y=`%IncMSE`)) +
  geom_segment( aes(x=reorder(Var.Names, `%IncMSE`), xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


#### Total consumption explanation ####


#### FNB TEST ####
test1 = read.csv("fnb-test.csv")
test1 = test1 %>%
  select(-e)

test2 = melt(test1, id.vars = c("session"))

test2 %>%
  group_by(session) %>%
  summarise(H = vegan::diversity(value, index = "shannon"),
         D = vegan::diversity(value, index = "simpson"))

test3 = test2 %>%
  group_by(session) %>%
  mutate(total = sum(value)) 

test3 = test3 %>%
  group_by(session) %>%
  rowwise() %>%
  mutate(prop = (value / total))

# fnb = food niche breadth
test3 = test3 %>%
  group_by(session) %>%
  summarise(
    sum.prop = sum(prop^2),
    Pi = (sum.prop^-1)-1, 
    N = n(),
    fnb = Pi/(N-1))

