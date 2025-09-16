# Filogenia
Se *[nada na biologia faz sentido exceto à luz da evolução](https://doi.org/10.2307/4444260)*, então é bastante justo que a discussão sobre como varia a forma seja feita sobre uma hipótese filogenética. Nas aulas [aula 7 e 8](Aulas%207%20e%208.pdf) falamos sobre isso, trabalhando conceitos como filomorfoespaço, sinal filogenético, convergência evolutiva, PGLS, taxas de evolução da forma e forma média. Usadas em conjunto, essas abordagens são complementares na interpretação e estudo da mudança da forma. A ideia aqui não é apenas visualizar como a filogenia se sobrepõe ao morfoespaço em si, mas quantificar e representar tendências evolutivas e lidar com a dependência filogenética em análises estatísticas. Em outras palavras: estamos quantificando a evolução morfológica.
Para essa prática, usaremos o mesmo conjunto da [aula 5](contour.md) (download [clicando aqui](Rhinella.TPS)), formado por landmarks `2D` posicionados sobre imagens de sapos do gênero *Rhinella* (retiradas da literatura). Na página *Curvas e contornos: semilandmarks e análise de Fourier* você pode encontrar mais detalhes.

## 1. Ocupação do filomorfoespaço e sinal filogenético
Antes de mais nada vamos instalar alguns pacotes adicionais que serão importantes para essa prática. O primeiro deles, `convevol`, calcula os índices C1, C2, C3 e C4 de [Stayton (2015)](https://doi.org/10.1111/evo.12729). Estes índices constituem métricas complementares de convergência morfológica, baseadas na tendência de aproximação entre terminais no morfoespaço em relação às posições de ancestrais comuns. Em outras palavras, são medidas do quanto diferentes linhagens tendem a ocupar uma região similar no morfoespaço, levando-se em conta como a forma de cada linhagem varia desde sua origem. Já o outro pacote, `phylocurve`, é usado para calcular o sinal filogenético através do índice *K<sub>mult</sub>*. O *K<sub>mult</sub>* é uma extensão multivariada do *K* de [Blomberg (2003)](https://doi.org/10.1111/j.0014-3820.2003.tb00285.x), proposta por [Adams (2014)](https://doi.org/10.1093/sysbio/syu030). Na sua versão original, o *K* é calculado como a razão entre a variância observada em uma característica quando ponderada pela estrutura filogenética e a variância calculada como se os valores na característica fossem independentes. Comentamos sobre isso na aula, você pode se aprofundar no tema acessando o site do [Prof. José Alexandre F. Diniz-Filho](https://dinizfilho.wixsite.com/dinizfilholab/sinal-filogenetico). O mais importante é ter em mente que se *K* = 1, então o traço varia exatamente como esperado pela filogenia; se *K* > 1, então o sinal filogenético é forte, porque a variãncia ponderada pelas relações de parentesco é maior do que a variância que ignora essas relações; e quando ocorre *K* < 1, ocorre o oposto e o sinal filogenético será baixo. Essa abordagem, no entanto, quantifica a relação com a filogenia apenas para traços contínuos e univariados. Já o *K<sub>mult</sub>* quantifica a intensidade do sinal filogenético em dados contínuos multivariados, como é o caso da posição ocupada por uma espécie no morfoespaço (i.e., descrita por valores específicos para cada PC). Se tiverem dificuldade na instalação do `phylocurve`, tentem acessar [a página do autor](https://ericgoolsby.github.io/software.html).
```{r instalar}
# Definir o diretório de trabalho
setwd("C:/caminho/para/pasta/desejada")

# Instalando pacotes novos
#install.packages("convevol")
#remotes::install_github("ericgoolsby/phylocurve") #se houver erro, instale o pacote devtools antes
```
Ótimo! Feito isso, podemos seguir com o carregamento normal dos pacotes e dados:
```{r carregar}
# Carregar pacotes necessários
library(geomorph)
library(tidyverse)
library(phytools)
library(convevol)
library(phylocurve)

# Carregar os dados
land.dt<-readland.tps("Rhinella.TPS", specID = "imageID", readcurves = TRUE)

# Verificando o número de dimensões do nosso tps
dim(land.dt)
```
A maneira mais simples de relacionar a ocupação do morfoespaço à filogenia é através da realização de um `filomorfoespaço`. Esse será o nosso primeiro objetivo no R e, para isso, vamos precisar: `1)` ler uma filogenia ultramétrica (e idealmente datada!) disponível na literatura (clique [aqui](TreePL-Rooted_Anura_bestTree.tre) para fazer o download da filogenia de [Portik et al. 2023](https://doi.org/10.1016/j.ympev.2023.107907), a mais recente para Anura); e `2)` cortar a filogenia de modo a manter apenas os terminais para os quais temos espécies representadas (no caso, as do gênero *Rhinella* que podem ser acessadas através de `dimnames(land.dt)[[3]]`).

```{r filogenia}
# Lendo a filogenia
tree<-read.tree("TreePL-Rooted_Anura_bestTree.tre")

# Se quisermos ver as espécies que coincidem entre os dois objetos
intersect(dimnames(land.dt)[[3]], tree$tip.label)
# Sequsermos ver as espécies que não coincidem
setdiff(dimnames(land.dt)[[3]], tree$tip.label)

# Atualizaremos o objeto tree para manter apenas as espécies desejadas na filogenia 
tree <- drop.tip(tree, setdiff(tree$tip.label, 
                         dimnames(land.dt)[[3]]))
```
Basicamente, estamos selecionando as espécies que não coincidem entre a filogenia e o `.tps` com `setdiff()` e removendo-as com `drop.tip()`. Você pode rodar `plotTree(tree)` para ver a árvore antes e depois desse processo (a diferença é bastante perceptível). Feito isso, ambos os objetos terão conteúdos referentes às mesmas espécies. No entanto, estão em ordens diferentes. No `R`, isso é um problema para quaisquer análises que envolvam o uso de filogenias. E como não podemos alterar a ordem das espécies na filogenia, vamos reordenar o nosso `.tps` para que esteja adequado:
```{r reordem}
# Ordenar o conjunto de dados para que coincida com a ordem das espécies na filogenia
land.dt <- land.dt[,,match(tree$tip.label, 
                           dimnames(land.dt)[[3]])]
```
Ótimo! Imagino que a essas alturas você já saiba como consultar a nova ordem das espécies no `.tps` (usando `dimnames(land.dt)[[3]]`). Por fim, podemos rodar a análise de Procrustes e gerar o morfoespaço.
```{r procrustes}
# Defina os sliders para cada curva
c1 <- define.sliders(12:42)
c2 <- define.sliders(41:52)
c3 <- define.sliders(51:82)
c4 <- define.sliders(81:101)

# Análise de Procrustes
gpa <- gpagen(land.dt, curves = rbind(c1, c2, c3, c4), 
              ProcD = FALSE)

# PCA para gerar o morfoespaço
pca <- gm.prcomp(gpa$coords)
```
Perfeito. Agora podemos gerar, por fim, um filomorfoespaço:
```{r filomorfoespaco1}
filomsp<-phylomorphospace(tree, X = pca$x[,1:2],
                 bty="n",node.size=c(0,1.2), label = "horizontal")
```
<p align="center">
<img src="filomorfspc1.png" alt="Fig1" width="630" height="490">
</p>

Legal, agora temos uma primeira representação do filomorfoespaço. Quando espécies filogeneticamente aparentadas ocupam uma região próxima no filomorfoespaço, como *R. ornata* e  *R. inopina* no nosso exemplo, podemos esperar uma tendência à manutenção da forma nas linhagens envolvidas. Se esse mesmo padrão estiver consolidado ao longo de toda a topologia, então esperamos encontrar um alto sinal filogenético (ou seja, *K<sub>mult</sub>* ≥ 1). 

```{r kmult}
# Estabelecendo um modelo evolutivo
evo.m <- evo.model(tree = tree,
                        Y = pca$x[,1:3],method = "Pairwise ML")

# Rodando o Kmult
Kmult<-K.mult(model = evo.m, nsim = 9999, plot = T)
```

## 2. Ecomorfologia e convergência
No nosso caso, o sinal filogenético foi relativamente baixo, o que significa que a similaridade morfológica não parece ser dependente da origem em comum entre linhagens. Por outro lado, ela pode ser resultado de seleção atuante na obtenção de recursos. Entre os anuros é comum que a morfologia reflita diretamente padrões ecológicos, representando um bom grupo para estudos ecomorfológicos. É como se determinada *forma* desempenhasse um papel primordial na *obtenção de recursos* em um determinado microambiente, de modo que exista uma tendência à maior similaridade entre espécies que ocupam aquele mesmo microambiente. Em outras palavras, pode haver uma convergência evolutiva atuante. Vamos testar algo nesse sentido, considerando que as espécies de *Rhinella* podem estar distribuídas em ambientes florestados ou de áreas secas. Para isso, criaremos um vetor indicando essa informação e extrairemos os dados do objeto `filomsp` para utilização no `ggplot()`:
```{r filomorfoespaco2}
# Vamos gerar um vetor atribuindo nossas categorias na mesma ordem em que as espécies estão dispostas no .tps
eco <- c("Mata", "Seca", "Mata", "Seca", "Seca", "Mata", 
         "Seca", "Seca", "Seca", "Mata", "Mata", "Mata", "Seca")
names(eco)<-tree$tip.label

# Extraindo do filomsp as informações de início e fim de cada segmento da filogenia para usar no ggplot 
filomsp.dt<-data.frame(
  xstart=filomsp$xx[filomsp$edge[,1]],
  ystart=filomsp$yy[filomsp$edge[,1]],
  xstop=filomsp$xx[filomsp$edge[,2]],
  ystop=filomsp$yy[filomsp$edge[,2]],
  nodestart=filomsp$edge[,1],
  nodestop=filomsp$edge[,2])

# Filomorfoespaço colorido por tipo de ambiente
ggplot()+
  geom_segment(data=filomsp.dt,aes(x=xstart,y=ystart,
                xend=xstop,yend=ystop), 
               linewidth = 0.75, colour="darkgray")+
  geom_point(data=pca$x[,1:2],
             aes(x=Comp1,y=Comp2, colour=eco),
             size=4, show.legend = T)+
  scale_color_manual(values = c("gold", "coral4"))+
  xlab("PC1")+
  ylab("PC2")+
  theme_bw()
```
<p align="center">
<img src="filomorfspc2.png" alt="Fig2" width="567" height="441">
</p>

Agora podemos testar se existe convergência morfológica entre as espécies em algum desses ambientes. Para isso, usaremos a função `convSig`, filtrando as linhagens para cada um dos grupos individualmente:

```{r convergencia}
# Testando a convergência morfológica entre as espécies que ocorrem em áreas florestadas
conv.Mata<-convSig(phy = tree, 
                   traits = as.matrix(pca$x[,1]), 
        focaltaxa = names(eco)[eco == "Mata"])

# Testando a convergência morfológica entre as espécies que ocorrem em áreas secas
conv.Seca<-convSig(phy = tree, 
                   traits = as.matrix(pca$x[,1]), 
        focaltaxa = names(eco)[eco == "Seca"])
```
Os resultados podem ser acessados através dos objetos `conv.Mata` e `conv.Seca`, e neles é possível observar o valor para cada um dos índices de Stayton (C1, C2, C3 e C4) e um `p-valor`. Cada um dos índices de Stayton quantifica um aspecto diferente da convergência evolutiva. Os índices de Stayton avaliam a redução da distância entre táxons ao longo do tempo, indicando a proporção dessa redução (C1), a magnitude absoluta da convergência morfológica (C2) e uma normalização desse valor em relação ao total de mudança evolutiva ao longo das duas linhagens (C3) ou em toda a filogenia (C4). Em outras palavras, os índices C1 e C2 se valem apenas de interpretações com base na distância de Procrustes para quantificar a convergência, enquanto C3 e C4 levam em consideração o tempo e os comprimentos de ramo. No nosso caso também não há muita convergência evolutiva, com uma redução absoluta (C1) de de aproximadamente 22% entre espécies de área florestada (`C1 = 0.2203`) e 26% entre espécies de áreas secas (`C1 = 0.2585`), sendo que estes valores não diferem do esperado ao acaso. Também podemos testar se a morfologia em cada ambiente tende a ser diferente. Vamos fazer isso junto a um teste alométrico:
```{r pgls}
# Rodando um PGLS
m1<-procD.pgls(phy=tree, coords ~ log(Csize)*eco, 
             data = gpa, iter = 9999) 
summary(m1)
```
Aqui, o que realizamos foi um PGLS (Phylogenetic Generalized Least Squares), uma regressão que “corrige” os dados para a história evolutiva compartilhada entre espécies, permitindo testar associações entre traços sem inflar falsos positivos devido à não-independência filogenética. E assim como nos casos anteriores, não constatamos efeitos atuantes sobre a forma.

## 3. Forma média e disparidade morfológica
Já sabemos que nesse conjunto de dados a forma não diverge tanto entre espécies de ambiente seco e florestado. Mas poderia divergir. E nessas circunstâncias uma abordagem que seria interessante seria identificar essas formas visualmente, de modo a comparar em quais aspectos elas diferem entre si. Para enxergar visualmente essa divergência, podemos calcular a *forma média* dos indivíduos. Aqui podemos usar o conceito de média aritmética como referência: uma medida de tendência central, representando o centro de um conjunto de dados. A forma média (ou *forma consenso*) segue a mesma ideia, e pode ser acessada através da função `mshape()`. A situação **mais comum** em que podemos precisar da forma média é quando temos múltiplos `.tps` (ou seja, `n` indivíduos) para a mesma espécie e desejamos fazer uma comparação *entre espécies* usando a filogenia. Essa comparação precisa ser feita descrevendo apenas uma forma por espécie. Então quando temos **múltiplas formas** é necessário calcular sua ***forma média*** para assegurar que o número de indivíduos seja o mesmo número de espécies na filogenia. Aqui será um pouco mais simples, bastando filtrar através da terceira dimensão do vetor `gpa$coords` os números dos espécimes em que `eco == "Mata"` e aqueles nos quais `eco == "Seca"`, para gerar dois novos objetos.
```{r mshape}
# Forma média
mean_shape_mata <- mshape(gpa$coords[ , , eco == "Mata"])
#plot(mean_shape_mata)
mean_shape_seca <- mshape(gpa$coords[ , , eco == "Seca"])
#plot(mean_shape_seca)

# Para visualizar as formas médias, vamos respeitar o padrão de coloração usado anteriormente. Para controlar as cores e os tamanhos dos pontos, usaremos
# a função gridPar()
par(mfrow = c(2,2))  #Para gerar um único plot
gp <- gridPar(grid.lwd = 0.25,
              grid.col = "lightgray",
              tar.pt.size = 1.1,
              tar.pt.bg=rep("gold",32)) #Usado para definir a cor das espécies de área florestada
plotRefToTarget(M1 = mean_shape_mata,
                M2 = pca$shapes$shapes.comp1$min, gridPars=gp) #Como é a forma no mínimo do PC1
plotRefToTarget(M1 = mean_shape_mata,
                M2 = pca$shapes$shapes.comp1$max, gridPars=gp) #Como é a forma no máximo do PC1

gp <- gridPar(grid.lwd = 0.25,
              grid.col = "lightgray",
              tar.pt.size = 1.1,
              tar.pt.bg=rep("chocolate",32)) #Usado para definir a cor das espécies de área seca
plotRefToTarget(M1 = mean_shape_seca,
                M2 = pca$shapes$shapes.comp1$min, gridPars=gp) #Como é a forma no mínimo do PC1
plotRefToTarget(M1 = mean_shape_seca,
                M2 = pca$shapes$shapes.comp1$max, gridPars=gp) #Como é a forma no máximo do PC1
``` 
<p align="center">
<img src="mshape.png" alt="Fig3">
</p>

Do mesmo modo que existe uma *medida de tendência central* para a forma, faz sentido pensar em como a forma *varia* dentro da amostra. Ou seja, deve haver alguma medida de *variãncia* da forma. E ela de fato existe, se chama *variância de Procrustes*. Vamos pensar num exemplo. Se pensarmos em um dado grupo distribuído no morfoespaço (indivíduos que formam uma espécie, indivíduos amostrados numa mesma localidade, diferentes espécies que compõem um gênero ou família, etc.), a `forma média` dá a ideia de ponto equidistante entre os indivíduos que compõem aquele grupo. A `variância de Procrustes` é a média das distâncias entre as formas dos indivíduos de um dado grupo e sua `forma média`. Em outras palavras, a ideia aqui é medir o quanto a nuvem de pontos que compõe esse mesmo grupo está mais ou menos "espalhada". Se imaginarmos que o morfoespaço é formado por `n` componentes principais (i.e. `PC1`, `PC2`, `PC3`, `PC4`...), então é fácil perceber que diversas medidas diferentes podem ser usadas para quantificar esse "espalhamento". A *área* ocupada por um grupo pode ser uma medida útil, quando utilizamos um morfoespaço formado por duas dimensões (`PC1` x `PC2`). O *volume* pode ser uma medida útil quando utilizamos mais dimensões que isso. Em conjunto, essas (`variância de Procrustes`, `área/volume ocupado por um grupo no morfoespaço`) e outras abordagens (e.g. diferença máxima entre as formas de um grupo no morfoespaço, ou `amplitude de Procrustes`) são medidas da ***disparidade morfológica***. No `R`, temos a função `morphol.disparity()` para verificar a variância de Procrustes diretamente.
```{r morphdisp1}
# Disparidade morfológica simples
m.disp<-morphol.disparity(gpa$coords ~ eco, 
                          groups = eco, transform = F)
summary(m.disp)
```
Veja que temos a definição *à priori* dos nossos grupos de interesse (o objeto `eco`), e através de uma fórmula indicamos o que esperamos: que a forma (acessada através do objeto `gpa$coords`) tenha variância distinta entre esses grupos. Aqui estamos usando uma fórmula simples, mas podemos aplicar um `PGLS` caso seja necessário corrigir a dependência filogenética entre os dados:
```{r morphdisp2}
# Disparidade morfológica com correção pelas relações de parentesco
m.disp.phy<-morphol.disparity(procD.pgls(gpa$coords ~ eco, phy = tree)
                          groups = eco, transform = F)
summary(m.disp.phy)
```
Em ambos os casos, a saída da função não só gera o valor da variância para cada um dos grupos, como também testa a diferença entre essas variâncias (conferindo um `p-valor` à análise).

## 4. Diversificação da forma no tempo evolutivo
Veja que, no fim das contas, a `forma média` e a `disparidade morfológica` são maneiras complementares de acessar a ocupação do morfoespaço, permitindo sua análise direta e quantitativa. Por outro lado, tínhamos tratado anteriormente o conceito de `filomorfoespaço`, em que a filogenia é usada para demonstrar relações entre as formas estudadas. Portanto, a *posição* que uma forma ocupa no morfoespaço pode ter significado biologicamente interpretável. E isso também é válido para formas *que são atualmente desconhecidas*. Os ramos da filogenia no filomorfoespaço, por exemplo, cruzam uma série de posições que podem ter sido ocupadas por formas extintas. Como já sabemos pela sistemática filogenética, filogenias são hipóteses, não uma realiadde absoluta. O mesmo vale aqui: é possível estabelecer hipóteses de evolução fenotípica examinando o morfoespaço. 
A primeira abordagem que usaremos nesse sentido é a reconstrução da forma ao longo da filogenia. Para isso, vamos reconstruir o valor do `PC1` (disponível em `pca$x[,1]`), já que é o melhor descritor da forma na nossa amostra, além de ser uma variável contínua. O primeiro passo é encontrar qual seria o valor esperado para os ancestras no `PC1`, o que é feito com a função `fastAnc()`. Aqui, a reconstrução será feita considerando um modelo de movimento browniano, para assumir o mínimo de pressupostos possível. Outros modelos podem ser adotados com funções específicas, e testes podem ser feitos para verificar qual modelo melhor se ajusta aos seus dados, mas isso merece uma página e prática próprias. 
```{r fastAnc}
PC1.rec<-c(pca$x[,1],fastAnc(tree,pca$x[,1]))
PC1.rec
```

Ao rodar `PC1.res` você deve ter observado que temos agora um valor para cada terminal e um valor para cada nó na árvore. Como a forma está sendo tratada como variável contínua, podemos extrapolar os valores esperados ao longo dos ramos utilizando a função `contMap`:
```{r contMapPC1}
# Mapeando a evolução da forma ao longo dos ramos: PC1
PC1.map <- contMap(tree, pca$x[,1], plot = FALSE)

# Visualizando o resultado da reconstrução
plot(PC1.map, legend = 0.7*max(nodeHeights(tree)), fsize = 0.7)

```
<p align="center">
<img src="reconst.png" alt="Fig4">
</p>

Perfeito Agora se você comparar a variação ao longo da filogenia mapeada vai perceber que estamos falando da mesma escala de valores no `PC1` no filomorfoespaço, de modo que pode interpretar como seria a forma esperada para cada nó. Também pode encontrar alguns padrões, como uma mudança importante na forma a partir do ancestral de ((*R. inopina*, *R. ornata*), *R. cerradensis*). De modo complementar, ao olhar as grades de distorção podemos inferir que essa mudança foi representada por uma tendência ao alongamento antero-posterior da glândula paratóide nessas espécies, acompanhada de lateralização da crista loreal e posicionamento mais anterior dos olhos. Obviamente, tudo isso também pode ser visualizado através do plot da grade de distorção do ancestral do grupo em relação à forma média. Para isso, temos que reconstruir o valor do `PC2`
```{r contMapPC2}
# Mapeando a evolução da forma ao longo dos ramos: PC2
PC2.rec <- c(pca$x[,2], fastAnc(tree, pca$x[,2]))
PC2.map <- contMap(tree, pca$x[,2], plot = FALSE)

# Encontrando o valor no PC1 e PC2 para o ancestral desejado
anc_PC1 <- PC1.rec[getMRCA(tree, c("Rhinella_inopina", 
                                   "Rhinella_ornata", 
                                   "Rhinella_cerradensis"))]

anc_PC2 <- PC2.rec[getMRCA(tree, c("Rhinella_inopina", 
                                   "Rhinella_ornata", 
                                   "Rhinella_cerradensis"))]
```
Os valores que estão "guardados" em `anc_PC1` e `anc_PC2` *são os mesmos* que a função `phylomorphospace()` gerou para o objeto `filomsp`, referenciando a posição do nó nº 16 (ou seja, o ancestral que estamos buscando; você pode confirmar isso rodando `filomsp$xx[16]` e `filomsp$yy[16]`). Isso mostra que todo o processo de elaboração do filomorfoespaço é matemático, utilizando a topologia da árvore e seus comprimentos de ramos na realização das estimativas. É por esse motivo que utilizar uma árvore datada é o correto nesse tipo de análise.


## 5. Inspecionando a forma nos 3 primeiros PCs
