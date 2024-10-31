# --- LOADING AND FORMAT -----------------------------------------------------
dat <- read.csv('data.csv')
head(dat)    # Una rapida occhiata al dataset
str(dat)     # Struttura del dataset

# Visualizzazione tabellare
kable(head(dat)) %>%
  kable_styling("striped", full_width = FALSE) %>%
  kable_paper() %>%
  save_kable(file = 'tableX.png', bs_theme = "flatly")

# --- PREPROCESSING ----------------------------------------------------------
sum(table(dat$Name))
colSums(is.na(dat))  # Controllo valori NA
dat$X <- NULL        # Rimuove colonna indici
dat$New_Price <- NULL
dat <- na.omit(dat)

# --- FORMATTAZIONE DELLA COLONNA "COMPANY" ----------------------------------
get_first_word <- function(text) {
  words <- strsplit(text, " ")[[1]]
  return(ifelse(length(words) > 0, words[1], NULL))
}
dat$company <- sapply(dat$Name, get_first_word)

# Separazione numerica per Mileage, Engine, Power
dat$Mileage <- c(unlist(sapply(dat$Mileage, get_first_word)), 0, 0)
dat$Engine <- unlist(sapply(dat$Engine, get_first_word))
dat$Power <- unlist(sapply(dat$Power, get_first_word))
dat <- na.omit(dat)

# Standardizza il tipo di "Owner"
dat[dat$Owner_Type == "Third", "Owner_Type"] <- "Third & Above"
dat[dat$Owner_Type == "Fourth & Above", "Owner_Type"] <- "Third & Above"

# --- ENCODING VARIABILI ORDINATE ---------------------------------------------
dat <- dat %>%
  mutate(
    Mileage = as.numeric(Mileage),
    Engine = as.numeric(Engine),
    Power = as.numeric(Power),
    Location = factor(Location),
    Fuel_Type = factor(Fuel_Type),
    Transmission = factor(Transmission),
    Owner_Type = factor(Owner_Type, levels = c("Third & Above", "Second", "First"), ordered = TRUE),
    Seats = factor(Seats),
    company = factor(company)
  )

# Filtraggio e rimozione dei valori indesiderati
dat <- na.omit(dat)
dat <- dat[dat$Mileage != 0, ]
dat <- dat[dat$Fuel_Type != "Electric", ]
dat$Fuel_Type <- factor(dat$Fuel_Type)
dat <- dat[-which.max(dat$Kilometers_Driven), ]
dat <- dat[-which.min(dat$Engine), ]

# --- RIORGANIZZAZIONE E VISUALIZZAZIONE DATI --------------------------------
num_dat <- dat[sapply(dat, is.numeric)]
chr_dat <- dat[sapply(dat, is.factor)]

# --- EDA --------------------------------------------------------------------
ggplot(dat, aes(x = Price)) + 
  geom_density(fill = "#E69F00", color = "black", alpha = 0.4) +
  geom_vline(xintercept = c(mean(dat$Price), median(dat$Price)), color = c("red", "blue"), linetype = "dotted", size = 0.5) +
  annotate("text", x = 25, y = 0.1, label = "Media Prezzo 9.6", color = "red") +
  annotate("text", x = 25, y = 0.13, label = "Mediana Prezzo 5.75", color = "blue") +
  labs(x = "Prezzo", y = "Densità", title = "Densità del Prezzo") +
  theme_minimal()


# --- MODELLI DI REGRESSIONE -------------------------------------------------
library(caret)
library(data.table)
library(mltools)
index <- createDataPartition(dat$Price, times = 1, p = 0.25)
test <- dat[index$Resample1, ]
train <- dat[-index$Resample1, ]

# Rimozione variabili non rilevanti
train$company <- NULL
test$company <- NULL

train_o <- one_hot(data.table(train))
test_o <- one_hot(data.table(test))

train_o$Owner_Type <- as.integer(ifelse(train$Owner_Type == "First", 1, ifelse(train$Owner_Type == "Secondo", 2, 3)))
test_o$Owner_Type <- as.integer(ifelse(test$Owner_Type == "First", 1, ifelse(test$Owner_Type == "Secondo", 2, 3)))

# Modello Gamma completo
gamma.fit.completo <- glm(Price ~ ., data = train_o, family = Gamma(link = "log"))
summary(gamma.fit.completo)

# Variabili significative per modello ridotto
vsign <- summary(gamma.fit.completo)$coeff[-1, 4] < 0.05
vsign <- names(vsign)[vsign == TRUE]
formulas <- as.formula(paste("Price ~", paste(vsign, collapse = "+")))

# Modello Gamma ridotto
gamma.fit.ridotto <- glm(formulas, data = train_o, family = Gamma(link = "log"))
summary(gamma.fit.ridotto)

# Test di likelihood ratio
anova(gamma.fit.ridotto, gamma.fit.completo, test = "LRT")

# Comparazione AIC e BIC per i modelli Gamma e Gaussian Inverso
AIC(gamma.fit.ridotto)
BIC(gamma.fit.ridotto)

# --- VISUALIZZAZIONE RISULTATI PREDITTIVI -----------------------------------
prediciton <- data.frame(
  P_gamma = exp(predict(gamma.fit.ridotto, test_o, type = "link")),
  P_igaussian = exp(predict(ig.fit.ridotto, test_o, type = "link")),
  test = test$Price
)

#MAE
mean(abs(P_gamma-test$Price))
mean(abs(P_igaussian-test$Price))

#RMSE
sqrt(mean((P_gamma-test$Price)^2))
sqrt(mean((P_igaussian-test$Price)^2))

#MSE
mean((P_gamma-test$Price)^2)
mean((P_igaussian-test$Price)^2)

# Plot confronto previsione Gamma vs Inverse Gaussian
p_gamma<-ggplot(prediciton, aes(x = test, y = P_gamma)) +
         geom_point(shape = 21, size = 1, alpha = 0.7, fill="blue") +  # Utilizza geom_point per un grafico a dispersione
         geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
         #scale_y_continuous(limits = c(0, 120))+  # Imposta i limiti dell'asse y
         labs(x = "Price_Test",y = "Prediction_gamma") +
         theme_minimal()

p_ig<- ggplot(prediciton, aes(x = test, y = P_igaussian)) +
        geom_point(shape = 21, size = 1, alpha = 0.7, fill="green") +  # Utilizza geom_point per un grafico a dispersione
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
        #scale_y_continuous(limits = c(0, 120))+  # Imposta i limiti dell'asse y
        labs(x = "Price_Test",y = "Prediction_igaussian") +
        theme_minimal()
grid.arrange(p_gamma,p_ig,nrow=1,ncol=2, top ="Confronto Test-Previsoni")


library(ggeffects)
ggpredict(gamma.fit.ridotto)
  

p<-ggpredict(gamma.fit.ridotto,terms = c("Kilometers_Driven","Power", "Transmission_Automatic"))
plot(p)


p<-ggpredict(gamma.fit.ridotto,terms = c("Power[50,100,150,200,300]","Owner_Type[1,2,3]","Year[2000,2008]"))
plot(p)



p<-ggpredict(gamma.fit.ridotto,terms = c("Kilometers_Driven","company_Fiat","Year[2010,2016]"))
plot(p)


p<-ggpredict(gamma.fit.ridotto,terms = c("Year"))
