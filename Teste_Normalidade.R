
#############################################################
#                     Geraração Nº Aleatórios               #
#############################################################

#############################################################
# Pacotes
library(DistributionTest)
library(moments)
library(nortest)
library(DescTools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
############################################################


############################################################
# Tamanho de Amostra
Resultados <- NULL
NN <- 10
#NN1 <- 20
#NN2 <- 50
#NN3 <- 100
############################################################


############################################################
escopo <- expand.grid(
  tamanho_amostra = seq(250, 300, 5),
  repeticao = seq(1, NN)
)

################################################################
Resultados <- purrr::map2_dfr(escopo$tamanho_amostra,escopo$repeticao, function(tamanho_amostra, ii){
  print(paste0(tamanho_amostra, "-", ii))
  #x <-rcauchy(tamanho_amostra, 0,1)
  #x <-rnorm(tamanho_amostra, 0,1)
  x <- rgamma(tamanho_amostra, shape = 10, rate = 1/10)
  tibble(
    Kolmogorov_Smirnov = ks.test(x, pnorm, mean(x), sd(x))$p.value,
    Jarque_Bera = JarqueBeraTest(x)$p.value,
    Anderson_Darling = AndersonDarlingTest(x, null = "pnorm", mean(x), sd(x))$p.value,
    Lilliefors = LillieTest(x)$p.value,
    Shapiro_Wilk = shapiro.test(x)$p.value,
    ZK = zk.test(x, 'norm')$p.value,
    ZC = zc.test(x, 'norm')$p.value,
    ZA = za.test(x, 'norm')$p.value,
    amostra = tamanho_amostra,
    tentativa = ii
  )
}
)
#######################################################################



##############################################################
# Comparação P-valor
Resultados |>
  tidyr::pivot_longer(
    names_to = "modelo",
    values_to = "estatistica",
    cols = c(Kolmogorov_Smirnov,
             Jarque_Bera,
             Anderson_Darling,
             Lilliefors,
             Shapiro_Wilk,
             #Cramer_Von_Mises,
             #agostino.test()
             ZK,
             ZC,
             ZA))|>
  group_by(amostra, modelo)|>
  summarise(
    minimo = min(estatistica),
    mediana = mean(estatistica),
    maximo = median(estatistica)+sd(estatistica))|>
  ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax =  maximo, 
             color = modelo))+
  geom_line()+
  geom_jitter()+
  geom_point()+
  ggtitle("Comparação do P-valo Médio")+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  labs(x = "Tamanho Amostral", 
       y = "P-Valor", 
       color = "Tipos de Testes")
#####################################################################
  

#####################################################################
# Poder do Teste
Resultados |>
  tidyr::pivot_longer(
    names_to = "modelo",
    values_to = "estatistica",
    cols = c(Kolmogorov_Smirnov,
             Jarque_Bera,
             Anderson_Darling,
             Lilliefors,
             Shapiro_Wilk,
             #Cramer_Von_Mises,
             #agostino.test()
             ZK,
             ZC,
             ZA))|>
  group_by(amostra, modelo)|>
  summarise(
    minimo = min(estatistica),
    mediana = mean(estatistica < .05),
    maximo = median(estatistica)+sd(estatistica))|>
  ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax =  maximo, color = modelo))+
  geom_line()+
  geom_jitter()+
  geom_point()+
  ggtitle("Comparação do Poder do Teste", subtitle = "Simulação")+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  labs(x = "Tamanho Amostral", y = "Poder do Teste", color = "Tipos de Testes")

  
###########################################################


DT::datatable(Resultados)





############################################################
# Distribuição Normal
Y <- rnorm(100000, 0,1)

# Distribuição Cauchy
Z <- rcauchy(100000)

# Distribuição Exponencial
W <- rexp(100000)

#Distribuição Uniforme
ZZ <- runif(1000)
#############################################################

#############################################################
# Visualização Gráfica
hist(Y)
hist(Z, breaks = 10)
hist(W)
hist(ZZ)
#############################################################

#############################################################
# Teste de Kolmogorov_Smirnov
ks.test(Y, "pnorm")
ks.test(Z, "pnorm")
ks.test(W, "pnorm")
ks.test(ZZ, "pnorm")
#############################################################

#############################################################
Amostra_de_p_valor_nulo <- numeric(length =  1000)
for(ii in seq(10,1000)){
  print(ii)
  Y_Amostra <- rnorm(10000, 0,1)
  Amostra_de_p_valor_nulo[ii] <- ks.test(Y_Amostra, "pnorm")$p.value 
}
#############################################################


#############################################################
# Distribuição Amostral do P-valor
#(Probability Integral Transform)

hist(Amostra_de_p_valor_nulo)

#############################################################


















