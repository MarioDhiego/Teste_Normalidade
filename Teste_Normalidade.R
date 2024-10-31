
################################################################################
#                     Geraração Nº Aleatórios                                  #
################################################################################

################################################################################
# Pacotes
library(DistributionTest)
library(irr)
library(moments)
library(nortest)
library(DescTools)
library(mvtnorm)
library(munsell)
library(fBasics)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
#library(kableExtra)
library(knitr)
library(tinytex)
library(Publish)
library(furrr)
library(rlecuyer)
################################################################################
################################################################################
# Tamanho de Amostra
Resultados <- NULL
NN <- 1000
#set.seed(1)
#set.seed(42)
################################################################################

################################################################################
# Geração de Dados
escopo <- expand.grid(
  tamanho_amostra = seq(25,30,5),
  repeticao = seq(1, NN)
)
################################################################################

################################################################################
# Simulaçãoe e cálculo dos testes estatísticos
Resultados <- purrr::map2_dfr(escopo$tamanho_amostra,escopo$repeticao, function(tamanho_amostra, ii){
  print(paste0(tamanho_amostra, "-", ii))
  #x <-rbeta(tamanho_amostra, 2,5)
  #x <-rcauchy(tamanho_amostra, 0,1)
  x <- rnorm(tamanho_amostra, 0,1)
  #x4 <- rgamma(tamanho_amostra, shape = 10, rate = 1/3)
  tibble(
    Kolmogorov_Smirnov = ks.test(x, pnorm, mean(x), sd(x))$p.value,
    Jarque_Bera = JarqueBeraTest(x)$p.value,
    Anderson_Darling = AndersonDarlingTest(x, null = "pnorm", mean(x), sd(x))$p.value,
    Lilliefors = LillieTest(x)$p.value,
    Shapiro_Wilk = shapiro.test(x)$p.value,
    Cramer_Von_Mises = cvm.test(x)$p.value,
    D_Agostino = agostino.test(x)$p.value,
    ZK = zk.test(x, 'norm')$p.value,
    ZC = zc.test(x, 'norm')$p.value,
    ZA = za.test(x, 'norm')$p.value,
    amostra = tamanho_amostra,
    tentativa = ii
  )
}
)
################################################################################



################################################################################
# Manipulaçãoe Visualização dos Resultados
# Comparação P-valor

Resultados %>%
  tidyr::pivot_longer(
    names_to = "modelo",
    values_to = "estatistica",
    cols = c(Kolmogorov_Smirnov,
             Jarque_Bera,
             Anderson_Darling,
             Lilliefors,
             Shapiro_Wilk,
             Cramer_Von_Mises,
             D_Agostino,
             ZK,
             ZC,
             ZA))%>%
  dplyr::group_by(amostra, modelo)%>%
  summarise(
    minimo = min(estatistica),
    mediana = mean(estatistica),
    maximo = median(estatistica)+ sd(estatistica))%>%
  ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax =  maximo,
             color = modelo))+
  geom_line()+
  geom_jitter()+
  geom_point()+
  #ggtitle("Comparação do P-valor Médio")+
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  labs(x = "Tamanho da Amostra(n)",
       y = "Erro Tipo I",
       color = "TESTES")
################################################################################


################################################################################
# Poder do Teste
Resultados %>%
  tidyr::pivot_longer(
    names_to = "modelo",
    values_to = "estatistica",
    cols = c(Kolmogorov_Smirnov,
             Jarque_Bera,
             Anderson_Darling,
             Lilliefors,
             Shapiro_Wilk,
             Cramer_Von_Mises,
             D_Agostino,
             ZK,
             ZC,
             ZA))%>%
  group_by(amostra, modelo)%>%
  summarise(
    minimo = min(estatistica),
    mediana = mean(estatistica < .05),
    maximo = median(estatistica)+sd(estatistica))%>%
  ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax =  maximo, color = modelo))+
  geom_line()+
  geom_jitter()+
  geom_point()+
  #ggtitle("Comparação do Poder do Teste",
   #       subtitle = "Simulação")+
  theme_bw()+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  labs(x = "Tamanho Amostral",
       y = "Poder do Teste",
       color = "Tipos de Testes")

# Gerar Tabela de Resultados
DT::datatable(Resultados)

mean(Resultados$Kolmogorov_Smirnov)
mean(Resultados$Anderson_Darling)
mean(Resultados$Shapiro_Wilk)

#------------------------------------------------------------------------------#
# Teste de Concordância
#Kappam.fleiss()

# Definir o limiar para a decisão binária (ajuste conforme necessário)
limiar <- 0.05

# Adicionar uma coluna de decisão binária para cada teste
Resultados_bin <- Resultados %>%
  mutate(
    Kolmogorov_Smirnov_bin = ifelse(Kolmogorov_Smirnov < limiar, 0, 1),
    Jarque_Bera_bin = ifelse(Jarque_Bera < limiar, 0, 1),
    Anderson_Darling_bin = ifelse(Anderson_Darling < limiar, 0, 1),
    Lilliefors_bin = ifelse(Lilliefors < limiar, 0, 1),
    Shapiro_Wilk_bin = ifelse(Shapiro_Wilk < limiar, 0, 1),
    Cramer_Von_Mises_bin = ifelse(Cramer_Von_Mises < limiar, 0, 1),
    D_Agostino_bin = ifelse(D_Agostino < limiar, 0, 1),
    ZK_bin = ifelse(ZK < limiar, 0, 1),
    ZC_bin = ifelse(ZC < limiar, 0, 1),
    ZA_bin = ifelse(ZA < limiar, 0, 1)
  )

# Calcular o Fleiss' Kappa para cada teste
kappa_results <- Resultados_bin %>%
  select(-c(amostra, tentativa)) %>%
  group_by(amostra) %>%
  summarise(across(everything(), ~ {
    kappa <- kappam.fleiss(as.matrix(.)[ ,1:8])  # Ajuste o número de avaliadores conforme necessário
    return(kappa$value)
  }
  )
  )




################################################################################




















################################################################################
# Distribuição Normal
Y <- rnorm(100000, 0,1)

# Distribuição Cauchy
Z <- rcauchy(100000)


# Distribuição gamma
x <- rgamma(1000, shape = 10, rate = 1/3)

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
hist(x)
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





#===========================================================


Normal_Classica <-rnorm(500,0,1)
Normal_Modificada <-rnorm(200,10,3)
t <-rt(500,35)
Exp<-rexp(500,2)
w<-rchisq(300,5)
f<-rf(500,10,12)

par(mfrow=c(3,2))
hist(Normal_Classica, col="lightblue4",border="white")
hist(Normal_Classica, col="lightblue4",border="white")
hist(t, col="lightblue4",border="white")
hist(Exp, col="lightblue4",border="white")
hist(w, col="lightblue4",border="white")
hist(f, col="lightblue4",border="white")


# Simulação
n = 5000
outc= sample(c("Head","Tail"), n, replace=T)
z = cumsum(outc=="Head")/seq(1,n)
plot(z, xlab="Flips", ylab="Frequency of Heads",type="l")
abline(h=0.5, col="grey")











