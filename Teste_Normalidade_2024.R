

# Carregar pacotes necessários
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
library(purrr)
library(knitr)
library(tinytex)
library(Publish)
library(furrr)
library(rlecuyer)
library(DT)
library(irr)






################################################################################
# Configurações
NN <- 2000  
tamanhos_amostra <- seq(25, 30, 5)  

################################################################################
# Criar escopo de simulação
escopo <- expand.grid(
  tamanho_amostra = tamanhos_amostra,
  repeticao = seq(1, NN)
)

################################################################################
# Simulação e cálculo dos testes estatísticos
Resultados <- purrr::map2_dfr(escopo$tamanho_amostra, escopo$repeticao, function(tamanho_amostra, ii) {
  print(paste0("Simulando: Tamanho Amostra = ", tamanho_amostra, ", Repetição = ", ii))
  #x <- rbeta(tamanho_amostra, 2,5)
  #x <- rnorm(tamanho_amostra, 0, 1)  # Gerar dados da distribuição normal
  x <- rcauchy(tamanho_amostra, location = 0, scale = 1)
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
})

################################################################################
# Calcular a TAXA MÉDIA de Acertos 
taxa_acertos <- Resultados %>%
  tidyr::pivot_longer(
    cols = c(Kolmogorov_Smirnov, 
             Jarque_Bera, 
             Anderson_Darling,
             Lilliefors,
             Shapiro_Wilk, 
             Cramer_Von_Mises, 
             D_Agostino, 
             ZK, 
             ZC, 
             ZA 
    ),
    names_to = "modelo",
    values_to = "estatistica") %>%
  dplyr::summarise(
    taxa_acertos = mean(estatistica > 0.05),  # Proporção de p-valores > 0.05
    .by = "modelo"                            # Agrupamento por "modelo"
  ) %>%
  dplyr::arrange(desc(taxa_acertos)) %>%
  dplyr::mutate(taxa_acertos = round(taxa_acertos, 3))

################################################################################
# Exibir a Tabela interativa para a Taxa Geral de Acertos
DT::datatable(
  taxa_acertos,
  options = list(pageLength = 10),
  caption = "Taxa Geral de Acertos por Teste Estatístico"
)

################################################################################
# Calcular a variância para cada teste
variancia_resultados <- Resultados %>%
  tidyr::pivot_longer(
    cols = c(Kolmogorov_Smirnov, 
             Jarque_Bera, 
             Anderson_Darling, 
             Lilliefors,
             Shapiro_Wilk, 
             Cramer_Von_Mises, 
             D_Agostino, 
             ZK, 
             ZC, 
             ZA
    ),
    names_to = "modelo",
    values_to = "estatistica"
  ) %>%
  dplyr::group_by(modelo) %>%
  dplyr::summarise(
    variancia = var(estatistica),  # Calculando a variância dos p-valores
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(variancia)) %>%
  dplyr::mutate(variancia = round(variancia, 3))


# Exibir a tabela interativa com a variância
DT::datatable(
  variancia_resultados,
  options = list(pageLength = 10),
  caption = "Variância dos P-valores por Teste Estatístico"
)

################################################################################
# Visualizar os Resultados em Gráficos

# Comparação do ERRO TIPO I (P-valor médio, mínimo e máximo)
Resultados %>%
  tidyr::pivot_longer(
    cols = c(Kolmogorov_Smirnov, 
             Jarque_Bera, 
             Anderson_Darling, 
             Lilliefors,
             Shapiro_Wilk, 
             Cramer_Von_Mises, 
             D_Agostino, 
             ZK, 
             ZC, 
             ZA
    ),
    names_to = "modelo",
    values_to = "estatistica"
  ) %>%
  dplyr::group_by(amostra, modelo) %>%
  dplyr::summarise(
    minimo = min(estatistica),
    mediana = mean(estatistica),
    maximo = median(estatistica) + sd(estatistica),
    .groups = "drop"
  ) %>%
  ggplot2::ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax = maximo, color = modelo)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Tamanho da Amostra (n)",
       y = "Erro Tipo I",
       color = "Testes")

################################################################################
# Visualizar o PODER DO TESTE
Resultados %>%
  tidyr::pivot_longer(
    cols = c(Kolmogorov_Smirnov, Jarque_Bera, Anderson_Darling, Lilliefors,
             Shapiro_Wilk, Cramer_Von_Mises, D_Agostino, ZK, ZC, ZA),
    names_to = "modelo",
    values_to = "estatistica"
  ) %>%
  dplyr::group_by(amostra, modelo) %>%
  dplyr::summarise(
    mediana = mean(estatistica < 0.05),  # Proporção de p-valores < 0.05
    minimo = min(estatistica),
    maximo = median(estatistica) + sd(estatistica),
    .groups = "drop"
  ) %>%
  ggplot2::ggplot(aes(x = amostra, y = mediana, ymin = minimo, ymax = maximo, color = modelo)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Tamanho da Amostra (n)",
       y = "Poder do Teste",
       color = "Testes")
################################################################################



