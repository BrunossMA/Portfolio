### Caricamento dei dati e ispezione iniziale #################################
dat <- read.csv('data.csv')
head(dat)   # Visualizza le prime righe del dataset per un'anteprima
str(dat)    # Struttura del dataset per verificare i tipi di variabili

### Preprocessing #############################################################
# Controllo per valori mancanti
colSums(is.na(dat))

# Conversione delle variabili categoriali in factor
dat <- dat %>%
  mutate(
    HeartDisease = factor(HeartDisease),
    Smoking = factor(Smoking),
    AlcoholDrinking = factor(AlcoholDrinking),
    Stroke = factor(Stroke),
    DiffWalking = factor(DiffWalking),
    Sex = factor(Sex),
    AgeCategory = factor(AgeCategory), 
    Race = factor(Race),
    Diabetic = factor(Diabetic),
    PhysicalActivity = factor(PhysicalActivity),
    GenHealth = factor(GenHealth),
    Asthma = factor(Asthma),
    KidneyDisease = factor(KidneyDisease),
    SkinCancer = factor(SkinCancer)
  )

# Rimozione di righe con eventuali NA residui
dat <- na.omit(dat)

### Analisi Esplorativa dei Dati (EDA) #######################################

# 1. Frequenza relativa dei soggetti con cardiopatia
df <- round(table(dat$HeartDisease) / nrow(dat), 2)
df <- data.frame("Status" = names(df), "Freq" = c(df[[1]], df[[2]]))

ggplot(df, aes(y = Freq, x = Status)) +
  geom_col(fill = c("red", "blue"), alpha = 0.4) +
  labs(
    x = "Cardiopatia", 
    y = "Frequenza relativa", 
    title = "Frequenza relativa dei soggetti con cardiopatia"
  ) +
  scale_x_discrete(labels = c("No", "Si")) +
  theme_minimal()

# 2. Distribuzione del BMI rispetto alla cardiopatia
ggplot(dat, aes(y = BMI, x = HeartDisease)) + 
  geom_boxplot(alpha = 0.4) +
  ylim(c(0, 100)) +
  labs(
    x = "Cardiopatia", 
    y = "Indice Massa Corporeo (BMI)", 
    title = "Box-plot BMI per Cardiopatia"
  ) +
  scale_x_discrete(labels = c("No", "Si")) +
  theme_minimal()

# 3. Attività fisica condizionata sulla presenza di cardiopatia
t <- table(dat$PhysicalActivity, dat$HeartDisease)
prop_tab <- prop.table(t, margin = 2)
prop_df <- as.data.frame(prop_tab)

ggplot(prop_df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Frequenza relativa di attività fisica condizionata su cardiopatia",
    x = "Attività Fisica", y = "Frequenza relativa"
  ) +
  scale_x_discrete(labels = c("No", "Si")) +
  scale_fill_discrete(labels = c("No", "Si")) +
  theme_minimal()

# 4. Frequenza di neoplasia della pelle condizionata su cardiopatia
t <- table(dat$SkinCancer, dat$HeartDisease)
prop_tab <- prop.table(t, margin = 2)
prop_df <- as.data.frame(prop_tab)

ggplot(prop_df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Frequenza relativa di neoplasia della pelle condizionata su cardiopatia",
    x = "Neoplasia", y = "Frequenza relativa"
  ) +
  scale_x_discrete(labels = c("No", "Si")) +
  theme_minimal()

### Correlazione e Selezione delle Variabili ##################################
# Creazione di un sottoinsieme di variabili numeriche
num_dat <- dat[sapply(dat[,-19], is.numeric)]
chr_dat <- dat[sapply(dat, is.factor)]
cor_matrix <- cor(num_dat)

# Visualizzazione della matrice di correlazione
corrplot(cor_matrix, method = "number", tl.col = "black", tl.srt = 45)

### Suddivisione Train-Test ##################################################
set.seed(1)
levels(dat$HeartDisease) <- c("0", "1")

# Creazione della partizione train-test
index <- createDataPartition(dat$HeartDisease, times = 1, p = 0.3)
test <- dat[index$Resample1, ]
train <- dat[-index$Resample1, ]

### Modelli di Regressione Logistica #########################################
# Modello di regressione logistica completo
log.comp <- glm(HeartDisease ~ ., data = train[,-19], family = binomial(link = "logit"))
summary(log.comp)

# Modello di regressione logistica ridotto
log.comp.rid <- glm(HeartDisease ~ ., data = train[, c(-19, -13)], family = binomial(link = "logit"))
summary(log.comp.rid)

# Confronto tra i due modelli logit
anova(log.comp, log.comp.rid, test = "LRT")

### Altri Modelli ############################################################
# Modello probit
prob <- glm(HeartDisease ~ ., data = train[,-19], family = binomial("probit"))
summary(prob)

# Modello cloglog
cll <- glm(HeartDisease ~ ., data = train[,-19], family = binomial("cloglog"))
summary(cll)

# Modello cloglog ridotto
cll.rid <- glm(HeartDisease ~ ., data = train[, c(-19, -13)], family = binomial("cloglog"))
summary(cll.rid)

# Confronto AIC tra modelli ridotti
AIC(log.comp.rid)
AIC(cll.rid)

# Calcolo delle previsioni e valutazione dei modelli
pred <- predict(log.comp.rid, newdata = test, type = "response")
pred <- factor(ifelse(pred > 0.5, "1", "0"), levels = c("0", "1"))
confusionMatrix(pred, test$HeartDisease)

### Ottimizzazione dei Pesi ##################################################
# Generazione di combinazioni di pesi per il modello pesato
step <- 0.01
n_weights <- 2
wmat <- generate_weights(step = step, nvar = n_weights)

# Inizializzazione del data frame per salvare i risultati
res <- data.frame(
  AIC = rep(0, nrow(wmat)), 
  Acc = rep(0, nrow(wmat)), 
  sensitivity = rep(0, nrow(wmat)), 
  specificity = rep(0, nrow(wmat)), 
  fraction_0 = wmat[,1], 
  fraction_1 = wmat[,2]
)

# Loop per valutare ogni combinazione di pesi
for (i in 1:nrow(wmat)) {
  train$V19 <- ifelse(train$HeartDisease == 0, wmat[i,1], wmat[i,2])
  log.pesato.rid <- glm(HeartDisease ~ ., data = train[, c(-19, -13)], weights = train$V19, family = binomial(link = "logit"))
  pred <- predict(log.pesato.rid, newdata = test, type = "response")
  pred <- factor(ifelse(pred > 0.5, "1", "0"), levels = c("0", "1"))
  
  a <- confusionMatrix(pred, test$HeartDisease)
  res[i, 1] <- AIC(log.pesato.rid)
  res[i, 2] <- a$overall[[1]]      # Accuratezza
  res[i, 3] <- a$byClass[[1]]      # Sensibilità
  res[i, 4] <- a$byClass[[2]]      # Specificità
}

print(res)
