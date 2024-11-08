# Tridimensionalidade
O primeiro "boom" de utilização da morfometria geométrica esteve atrelado à utilização de fotografias e esquemas para levantamento de dados em duas dimensões. Programas como `TpsUtil` e `TpsDig` tiveram um papel importante nesse processo, como já comentamos anteriormente. Mas até que ponto a representação da forma a partir de dados bidimensionais é uma boa tradução da realidade?  
A "perda da tridimensionalidade" tem sido alvo de debate teórico por décadas (veja [Adams et al., 2004](https://doi.org/10.1080/11250000409356545) e [Cardini, 2014](https://doi.org/10.4404/hystrix-25.2-10993) para detalhes a respeito), e é concomitante ao desenvolvimento da própria morfometria geométrica em si. Felizmente, tem existido uma ampliação de acesso a métodos que possibilitam a obtenção de dados 3D. Aqui podemos citar principalmente a *fotogrametria* e a *tomografia computadorizada*.  
Comentamos um pouco mais sobre eles na [aula 6](Aula%206.pdf). A prática em si foi focada no uso de dados tomografados, fizemos a reconstrução de uma tomografia em formato NIfTI (`.nii`) utilizando o software livre brasileiro [InVesalius](https://www.gov.br/cti/pt-br/acesso-a-informacao/acoes-e-programas/invesalius). Esse processo nos gera os modelos tridimensionais de superfície (as chamadas *meshs*). A seguir, utilizamos o IDAV Landmark Editor para posicionar landmarks e semilandmarks sobre algumas meshs de crânios de anuros (seguindo o tutorial disponível [aqui](Aula%206%20-%20Tutorial.pdf)). As *meshs* utilizadas foram obtidas a partir do [MorphoSource](https://morphosource.org/), base de dados bastante útil para nós *morfólogos*.   
Para dados tridimensionais, existe a possibilidade de aplicar os semilandmarks e um contexto mais amplo: ao invés de limitar-los à descrição de *curvas* também podemos utilizá-los para descrever ***superfícies***. Além disso, entre as *meshs* utilizadas está um fóssil, Electrorana limoae, anuro cretácico encontrado preservado no ambar em Myanmar (clique [aqui](https://doi.org/10.1038/s41598-018-26848-w) para a descrição, publicada na Nature). Isso nos dá a oportunidade de trabalharmos duas abordagens importantes: a **retrodeformação** e a estimativa de landmarks em casos de **missing data**. O conjunto de dados as coordenadas para todos os landmarks pode ser baixado [aqui](Aula_6-dta.zip).  

## 1. Semilandmarks de superfície
Vamos começar carregando os dados (que agora estão no formato `.dta`) e chamando os pacotes necessários. Lembre-se de instalá-los caso não os tenha.
```{r data}
# Definir o diretório de trabalho
setwd("C:/caminho/para/pasta/desejada")

# Carregar pacotes necessários
library(geomorph)
library(Morpho)
library(tidyverse)
library(ggpubr)

# Listando os arquivos .dta em um objeto...
lista<-list.files(pattern = ".dta")
# (...) e usando esse objeto para leitura
land.dt<-readmulti.nts(lista)
```
Se rodarmos `dimnames(land.dt)[[3]]`, vemos o nome das espécies que fazem parte do nosso dataset. Perceba que `Electrorana_limoae` ocupa a *terceira* posição. Essa informação será útil mais adiante. Por agora perceba que a complexidade das coisas tende a aumentar um pouco. De modo similar aos semilandmarks de curvas que retratamos na aula passada, semilandmarks de superfície também devem ser submetidos a `sliding`. Portanto, precisamos conseguir definir quais são os landmarks, semilandmarks que descrevem curvas e que descrevem superfícies. Em outras palavras, precisamos conhecer bem nosso conjunto de dados.  
Como foi comentado em sala de aula, quando rodamos `dim(land.dt)` para encontrar o número de dimensões do nosso objeto, obtemos `[1] 343   3   6`. O primeiro valor, `343`, é justamente referente ao número de marcos anatômicos utilizado. Sabemos que os primeiros 21 deles são landmarks padrão, fixos:  

```{r landmarks}
####Definindo landmarks e semis
#De 1 a 21, temos landmarks. Daí em diante, é tudo semi 
gpa<-gpagen(land.dt[1:21,,])
plotAllSpecimens(gpa$coords)
```

<p align="center">
<a href="https://752ynx-diego-almeida0silva.shinyapps.io/Aula6_Ldkfix/">
  <img src="Aula6_ldkfix.png" alt="Gráfico Interativo" title="Clique na figura" width="600" height="450">
</a>
</p>

Daí por diante, todos os demais são semilandmarks. Os primeiros correspondem às curvas, porque essa foi a ordem seguida durante o posicionamento no IDAV Landmark Editor. Se temos duas curvas de 30 pontos seguidas de outras quatro formadas por 20, então os próximos 140 pontos representam semilandmarks de curvas.

```{r landmarks}
####As nossas curvas são de 30+30+20+20+20+20 = 140ldk. 
#De 22 a 162, portanto, estamos descrevendo curvas
gpa<-gpagen(land.dt[22:162,,])
#gpa<-gpagen(land.dt)
plotAllSpecimens(gpa$coords)
```
<p align="center">
<a href="https://752ynx-diego-almeida0silva.shinyapps.io/Aula6_Curvas/">
  <img src="Aula6_curvas.png" alt="Gráfico Interativo" title="Clique na figura" width="600" height="450">
</a>
</p>

Já a nossa superfície é descrita por 200 pontos, começando a partir do 143 e indo até o último presente. Para tornar o script mais genérico, vamos definir que estamos selecionando os seilandmarks do objeto `land.dt` entre `163:dim(land.dt)[1]` (ou seja, partindo de 163 até o último ponto disponível).

```{r landmarks}
####As nossa superfície vai daí em diante: [143,dim(land.dt)[1]] = 200 ldk
gpa<-gpagen(land.dt[163:dim(land.dt)[1],,])
#gpa<-gpagen(land.dt)
plotAllSpecimens(gpa$coords)
```

<p align="center">
<a href="https://752ynx-diego-almeida0silva.shinyapps.io/Aula6_Superf/">
  <img src="Aula6_superficie.png" alt="Gráfico Interativo" title="Clique na figura" width="600" height="450">
</a>
</p>
