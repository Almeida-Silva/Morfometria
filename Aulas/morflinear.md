# Morfometria Linear 
Começamos as aulas [1](Aula%201.pdf) e [2](Aula%202.pdf) tratando de morfometria linear. A ideia aqui é contextualizar seu uso ao longo do tempo e discutir sobre as abordagens estatísticas mais adequadas para esse tipo de dado. Já a atividade prática envolveu o uso do programa [ImageJ](https://imagej.net/software/fiji/downloads) para a obtenção de medidas lineares a partir de [fotos de girinos disponíveis na literatura](Fotos_Aula1.zip). Acontece que a família Hylodidae, aqui usada como exemplo, possui dois gêneros (_Megaelosia_ e _Phantasmarana_) cujas larvas podem chegar a ± 12cm, e outros dois (_Crossodactylus_ e _Hylodes_) que apresentam girinos de ± 4cm. Será que as medidas lineares são suficientes pra sugerir que esses bichos diferem entre si, ou será que toda a diferença encontrada reside na discrepância do tamanho?
Para testar isso, usaremos os dados obtidos pelos alunos e salvos na planilha "[Dados_aula1.xlsx](Dados_aula1.xlsx)". 

## Resíduos das medidas lineares pelo tamanho
Uma abordagem bastante aceita envolve a remoção do efeito do tamanho através de uma propriedade estatística: os resíduos de um modelo. Vamos começar indicando nossa pasta de trabalho, carregando os pacotes e importando os dados diretamente de uma planilha do Excel. Se você abrir a planilha no computador, repare que as medidas lineares foram tomadas entre as colunas C e O, enquanto A e B são as colunas que designam a que espécie e gênero cada indivíduo pertence. 

```{r data}
# Definir o diretório de trabalho
setwd("C:/caminho/para/pasta/desejada")

# Carregar pacotes necessários
library(readxl)
library(missMDA)
library(nlme)
library(tidyverse)
library(MASS)

# Carregar os metadados
metadados <- read_xlsx("Dados_aula1.xlsx", sheet = 'Aula1', range = cell_cols("A:B"))

# Carregar os dados numéricos
dados <- read_xlsx("Dados_aula1.xlsx", sheet = 'Aula1', range = cell_cols("C:O"), col_types = 'numeric')
```

**Estimar os NA**  
Um problema comum na morfometria linear é a impossibilidade de tomar uma ou mais medidas de algum(ns) indivíduo(s). Existem algumas maneiras de lidar com esse tipo de questão, como remover a variável ou o espécime (quando isso não parecer impactar de maneira importante no conjunto de dados) ou buscar um método de inferência para indicar qual seria o valor esperado. Aqui aplicaremos essa segunda abordagem. Para isso, seguiremos a proposta de [Josse & Husson 2016](https://doi.org/10.18637/jss.v070.i01), um método que se baseia na análise de componentes principais. 

```{r missMDA, echo=FALSE, eval=FALSE}
# Estimar o número ótimo de componentes principais para imputar os dados faltantes
nb <- estim_ncpPCA(dados, ncp.min = 0, ncp.max = 5, method.cv = "Kfold", nbsim = 1000)
nb # Número ótimo de componentes principais

# Imputar os valores faltantes
res.impute <- imputePCA(dados, ncp = nb$ncp)
dados <- round(res.impute$completeObs, 1)
rownames(dados) <- metadados$spp
```

**Transformação por log**  
Uma vez que o conjunto de dados esteja completo, é necessário transformá-lo de modo que os dados estejam em uma mesma escala. Para isso, aqui transformaremos os dados através de seu logaritmo. Entre outras vantagens dessa abordagem, ela reduz o efeito de valores extremos na amostra e normaliza sua distribuição. É válido comentar que a forma como é usada a função `log()` abaixo converte os dados originais para seu logaritmo natural (neperiano), dado que não indicamos um valor para a base. Consulte o help da função para mais detalhes.

```{r data_v2, echo=FALSE, eval=TRUE}
# Transformar os dados com logaritmo natural
dados <- log(dados)
```

**Modelo de regressão linear (RL)**  
Dado que queremos remover o efeito do tamanho para analisar se continua existindo diferenciação entre as medidas morfológicas, vamos empregar aqui uma regressão linear múltipla. A ideia aqui se aproveita de uma propriedade estatística. Quando criamos um modelo para testar se uma variável *y* na natureza pode ser explicada por outra *x*, sabemos que essa relação não será perfeita: existe uma certa variação nos valores de *y* que está além do que *x* pode explicar. Esse distanciamento entre o que seria esperado em um modelo perfeito e o que de fato observamos nos dados recebe o nome de *resíduo*. Então se gerarmos um modelo em que testamos se todas as medidas são explicadas pelo tamanho, os resíduos desse modelo estarão minimamente livres do efeito do tamanho em si.  

```{r dividir, echo=FALSE, eval=TRUE}
# Definir variáveis de medidas e tamanho
medidas <- as.matrix(dados[, 2:ncol(dados)])
size <- dados[, 1]

# Ajustar um modelo de regressão
regr <- gls(medidas ~ size)
summary(regr)
```

**Análise de componentes principais (PCA)**  
Uma PCA é uma análise que visa reduzir a dimensionalidade dos dados através da maximização da explicação da variância, o que é feito encontrando eixos que representam combinações lineares das variáveis originais e que são perpendiculares entre si. [Esse vídeo](https://youtu.be/FgakZw6K1QQ) pode ser bastante esclarecedor para quem quiser entender melhor como funciona essa análise

```{r pca, echo=FALSE}
# PCA com os resíduos da regressão
pca.res <- prcomp(regr$residuals)
summary(pca.res)

# Preparar os dados para o gráfico
df.PCAres <- data.frame(spp = metadados$spp, 
                        gen = metadados$gen, 
                        PC1 = pca.res$x[, 1], 
                        PC2 = pca.res$x[, 2], 
                        PC3 = pca.res$x[, 3])

# Gráfico da PCA (PC1 vs PC2)
ggplot(df.PCAres, aes(PC1, PC2)) +
  geom_point(aes(fill = gen), color = "black", size = 4.5, pch = 21, show.legend = TRUE) +
  xlab("PC1") +
  ylab("PC2") +
  theme_bw()
```

## MANOVA

Existe de fato diferência entre os gêneros?

```{r teste, echo=FALSE}
# Teste MANOVA
teste <- manova(medidas ~ df.PCAres$gen)
summary(teste)
```


## Posthoc

Usando uma LDA para como se separam os grupos

```{r posthoc, echo=FALSE}
# Análise discriminante linear
lda.phc <- lda(df.PCAres$gen ~ medidas, CV = FALSE)

# Preparar os dados da LDA para o gráfico
lda_df <- data.frame(
  gen = df.PCAres[, "gen"],
  lda = predict(lda.phc)$x
)
lda_df

# Gráfico da LDA
ggplot(lda_df) +
  geom_point(aes(x = lda.LD1, y = lda.LD2, fill = gen), 
             color = "black", pch = 21, size = 4) +
  labs(x = "LD1", y = "LD2") +
  theme_bw()
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
