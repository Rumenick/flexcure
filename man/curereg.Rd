\encoding{latin9}
\name{curereg}
\alias{curereg}
\title{
  Regression for a flexible parametric survival
  models with cure fraction
  %Regressão para modelos flexíveis paramétricos de sobrevivência com fração de cura
}
\description{
  Ajusta modelos de regressão para distribuições
  flexíveis paramétricas de sobrevivência com
  fração de cura. Pode-se estimar, simultaneamente,
  o efeito das covariáveis sobre o tempo até a
  ocorrência do evento e fração de curados.
%Fit a flexible parametric survival regression models with cure fraction.
}


\usage{
curereg(formula, cureformula = ~1, data, weights, timedist = "moeev",
        ncausedist = "poisson", subset, na.action,
        inits, method = "SANN", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula }{
    expressão fórmula na sintaxe convencional do \code{R}, como para outros modelos de regressão
    . A resposta é um objeto \code{survival}, como retornado pela função \code{\link{Surv}} e as
    covariáveis são especificadas ao lado direito  de \code{~}, por exemplo:

    \code{formula = Surv(time, status) ~ age + sex}

    São suportados dados censurados à direita, à esquerda e intervalar. Para mais detalhes,
    consulte a documentação de \code{\link{Surv}} e \code{\link{flexsurvreg}}.
    Se as covariáveis não são utilizadas, então especifique \code{1} ao lado direito de \code{~}
    , por exemplo \code{formula = Surv(time, status) ~ 1}.
  }
  \item{cureformula }{
    expressão fórmula sem resposta, que é usada para especificar o efeito das covariáveis sobre
    a fração de cura, por exemplo:

    \code{cureformula = ~ age + sex}

    se as covariáveis não são utilizadas, então especifique \code{1} ao lado direito de \code{~}
    , por exemplo \code{cureformula = ~ 1} (default). Se \code{cureformula = } ou \code{cureform
    ula = NULL}, será ajustado o modelo sem fração de cura especificado em \code{timedist}
    usando diretamente a função \code{\link{flexsurvreg}}, exceto a distribuição valor extremo
    na forma estendida de Marshall-Olkin ("moeev"), pois não faz parte das distribuições
    disponíveis através desta função.
  }
  \item{data }{
    cojunto de dados (\code{data.frame}) ou lista (\code{list}), em que se encontram as
    variáveis definidas na \code{formula} e \code{cureformula}. Se \code{data} não for informado
    , as variáveis devem estar no ambiente de trabalho \code{(.GlobalEnv)}.
  }
  \item{weights }{
    Vetor opcional de pesos a serem utilizados no ajuste.
  }
  \item{timedist }{
    distribuição associada ao tempo de sobrevivência dos elementos não curados. Podem ser
    utilizadas as seguintes distribuições: \code{"exp"} (exponencial),
    \code{"weibull"} (Weibull
    ), \code{"ev"} (valor extremo), \code{"gamma"} (Gama),\code{"lnorm"} (Log-normal),
    \code{"llogist"} (Log-logística), \code{"moee"} (valor extremo (or exponencial) na forma estendida de Marshall-Olkin),
    \code{"moeev"} (valor extremo (or Weibull) na forma estendida de Marshall-Olkin), \code{"gengamma"}
    (Gama generalizada estendida), \code{"genf"} (F generalizada), sendo o padrão
    \code{timedist = "moeev"}.

    A parametrização usada para as distribuições exponencial, Weibull, log-normal e
    log-logística é diferende da utilizada em \code{survreg}. Para os modelos exponencial,
    Weibull e log-normal usa-se a parametrização definidas, respectivamente, pelas funções:
    \code{\link{dexp}}, \code{\link{dweibull}} e \code{\link{dlnorm}} do pacote \pkg{base}. Já
    para a distribuição log-logística é utlizada a parametrização definida pela função
    \code{\link{dllogist}} do pacote \pkg{flexsurv}, assim como também para os modelos gama
    generalizado estendido (\code{\link{dgengamma}}) e F generalizado (\code{\link{dgenf}}). O
    modelo valor extremo na forma extendida de Marshall-Olkin é definido em
    \code{\link{dmoeev}} no pacote \pkg{flexsurv.cure}.
  }
  \item{ncausedist }{
    distribuição associada ao número de causa que competem para ocorrência do evento, ou seja,
    se \code{ncausedist = "bernoulli"}, tem-se o modelo de mistura padrão e se
    \code{ncausedist = "poisson"} (default), o modelo de tempo de promoção.
  }
  \item{inits }{
    vetor númerico informando valores iniciais para os parâmetros desconhecidos, porém opcional.
    Estes devem ser informado na seguinte ordem, primeiro os valores iniciais do termo relaciona
    do com a fração de cura, em segundo os valores inicias do termo relacionado com o tempo de
    sobrevivência e, por fim, os valores iniciais para os demais parâmetros associados a
    \code{timedist}. Por exemplo, deseja-se ajustar o modelo de tempo de promoção gama generaliz
    ado com duas covariáveis no termo associado ao tempo, então

    \code{inits = c(intercetp_cure = .1, intercetp_time = mean(log(t)),
    cov1 = 0, cov2 = 0, sigma = sd(log(t)), Q = 0)}

    sendo os 'names' utilizados apenas para efeito ilustrativo.

    Se não for especificado, os valores iniciais padrão são escolhidos a partir de um simples
    resumo dos tempos de sobrevivência ou do indicador de censura, por exemplo, a média é muitas     vezes usada para inicializar parâmetros de posição e a proporção de censura para inicializar
    o parâmetro associado a fração de cura. Veja 'detalhes' e os objetos curereg.dists ou
    flexsurv.dists na fonte para os métodos exatos utilizados. Além disso, caso a superfície de
    log-verossimilhança seja irregular, é aconselhável executar a otimização a partir de
    diferentes valores iniciais, para assegurar a convergência para o verdadeiro máximo global.
  }
  \item{subset }{
    expressão informando que apenas um subconjunto das linhas dos dados deve ser utilizado no
    ajuste.
  }
  \item{na.action }{
    função de filtro para dados faltantes, em que é aplicada após \code{subset}, se for o caso.
    Por padão é usado "na.omit".
  }
  \item{method }{
    método de otimização a ser utilizado na função \code{optim} (rotina de otimização padrão do
    \code{R}). Por padrão \code{method = "BFGS"}, porém podem ser especificados também
    "Nelder-Mead", "CG", "L-BFGS-B" e "SANN". Para detalhes sobre os métodos consulte
    \code{\link{optim}}.
  }
  \item{\dots }{
    outros argumentos opicionais que podem ser passados para as funções \code{\link{optim}} e
    \code{\link{flexsurvreg}}.

    no caso da função \code{optim}, se a otimização não converge, considere normalizar o
    problema utilizando, por exemplo, control = list(fnscale = value), sendo \code{value} um
    número da ordem de magnitude da função de log-verossimilhaça. Além disso, quando a
    otimização não converge, a rotina libera uma messagem que a Hessiana não é positiva
    definida. Nesta situação, é sugerido que o usuário aperte os critérios de tolerância para a
    convergência. De modo geral, se a otimização leva muito tempo, os passos intermediários
    podem ser impressos usando o argumento \code{trace} na lista de controle. Veja
    \code{\link{optim}} para mais detalhes.

    no caso da função \code{flesurvreg}, quando deseja-se fixar algum parâmetro em um modelo a
    ser estimado, pode-se usar o argumento \code{fixedpars}. Em que este é um vetor de índices
    dos parâmetros, cujos valores serão fixados em seus valores iniciais durante a otimização.
    Os índices são ordenados como em \code{inits}. Por exemplo, em um modelo de mistura padrão
    gama generalizado com duas covariáveis no termo relacionado com fração de cura, para fixar o
    terceiro de quatro parâmetros 'base' deste modelo (parâmetro \code{Q}, consulte a ajuda de
    \code{\link{gengammasm}}) e o efeito da segunda covariável, especifique fixedpars = c(3, 6).
  }
}
\details{
  Note que os argumentos \code{ncausedist} e \code{timedist} define qual modelo será ajustado,
  ou seja, se \code{ncausedist = "poisson"} e \code{timedist = "genf"}, tem-se o modelo de
  tempo de promoção F generalizado. A função de densidade de probabilidade imprópria deste
  modelo está definida em \code{\link{dgenfpt}}, que pode ser acessada pelo usuário. Se
  \code{ncausedist = "bernoulli"} e \code{timedist = "genf"}, o modelo de mistura padrão é
  ajustado. Da mesma forma, o usuário pode acessar a função de densidade de probabilidade
  imprópria em \code{\link{dgenfms}}. Por definição, tem-se \code{dnamems} e \code{pnamems},
  para os modelos de mistura padrão, sendo 'name' o nome da distribuição da forma especificada
  em \code{timedist} e o sufixo \code{ms} que significa 'standart mixture'. Assim como, para os
  modelos de tempo de promoção, a única diferença é que o sufixo é o \code{pt}, que significa
  'promotion time'.

  As covariáveis são inseridas no termo relacionado ao tempo dos elementos não curados,
  considerando um modelo de tempo de falha acelerado (Lawless, 2003).  Já para o termo relaciona
  do com a fração de cura, se \code{ncausedist = "bernoulli"} as covariáveis são inseridas de
  maneira similar a um modelo de regressão logístico com função de ligação logística
  (\code{family = binomial(link = "logit")} em \code{glm}). Quando \code{ncausedist = "poisson"}
  as covariáveis são inseridas de maneira similar a um modelo regressão poisson com função de
  ligação logarítmica (\code{family = poisson(link = "log")} em \code{glm}). Para mais detalhes
  consulte as referências.
}

\value{
  Lista de classe "curereg" contendo informações sobre o modelo ajustado, que pode ser estudada
  através das funções: \code{print}, \code{summary}, \code{\link{coef.curereg}},
  \code{\link{curefraction}}, \code{\link{plot.curereg}} e \code{\link{lines.curereg}}.
  Componentes de interesse para os usuários podem incluir:

  \item{call }{
    cópia da função executada a ser usada pós-processamento.
  }
  \item{
    coefficients }{lista com as estimativas de máxima verossimilhança dos parâmetros
    (ver detalhes).
  }
  \item{std.error }{
    lista com as estimativas de erros padrão dos parâmetros (ver detalhes).
  }
  \item{vcov }{
    matriz de variância e covariância dos parâmetros, que pode ser obtida com
    \code{\link{vcov}} (ver exemplo).
  }
  \item{loglik }{
    log-verossimilhaça.
  }
  \item{AIC }{
    critério de informação Akaike (- 2 * log-verossimilhaça + 2 * número de parâmetros
    estimados).
  }
}
\references{
  Jackson, Christopher H. "flexsurv: a platform for parametric survival modelling in R."

  Maller, R. A., & Zhou, X. (1996). Survival analysis with long-term survivors.
  New York: Wiley.

  Lawless, J. F. (2011). Statistical models and methods for lifetime data (Vol. 362).
  John Wiley & Sons.

  Ortega, E. M., Cancho, V. G., & Paula, G. A. (2009). Generalized log-gamma regression models
  with cure fraction. Lifetime Data Analysis, 15(1), 79-106.

  Peng, Y., Dear, K. B., & Denham, J. W. (1998). A generalized F mixture model for cure rate
  estimation. Statistics in medicine, 17(8), 813-830.
}

\author{
  Rumenick Pereira da Silva <rumenickbf@hotmail.com>
}

\note{
%%  ~~further notes~~
}

\seealso{
  \code{\link{confint.curereg}} (estima intervalo de confiança para os parâmetros),
  \code{\link{predict.curereg}} (estima a fração de cura), \code{\link{plot.curereg}} (retorna
  gráfico da sobrevivência estimatida por Kaplan-Meier sobreposta a estimada pelo modelo
  ajustado) e \code{\link{lines.curereg}} (desenha sobre \code{plot.curereg} a sobrevivência
  estimada para outro modelo ajustado, útil para realizar comparações).
}
\examples{
## fit Extended generalized gamma standart mixture model

if(!any("smcure" == rownames(installed.packages()))){
install.packages("smcure", repos = "http://cran.us.r-project.org")}

require(smcure) # for data 'e1684' e 'bmt'
data(e1684)
fitgg <- curereg(Surv(FAILTIME, FAILCENS) ~ TRT + SEX + AGE, cureformula = ~ TRT + SEX + AGE,
                 data = e1684, timedist = "gengamma", ncausedist = "bernoulli")

# Output of 'curereg' object
fitgg
# Extract Model Coefficients
coef(fitgg)
# Extract model coefficients:
# Terms: failure time distribution model
coef(fitgg, terms = "time")
# Terms: cure probability model
coef(fitgg, terms = "cure")
# Object summaries
summary(fitgg)
# Information criterion
AIC(fitgg) # Akaike information criterion
BIC(fitgg) # Bayesian information criterion
# Extract Log-Likelihood
logLik(fitgg)
# Calculate Variance-Covariance Matrix for a Fitted Model Object
vcov(fitgg)

## fit Generalized F standart mixture model
fitgf <- curereg(Surv(FAILTIME, FAILCENS) ~ TRT + SEX + AGE, cureformula = ~ TRT + SEX + AGE,
                 data = e1684, timedist = "genf", ncausedist = "bernoulli")

# Output of 'curereg' object
fitgf
# Extract Model Coefficients
coef(fitgf)
# Extract model coefficients:
# Terms: failure time distribution model
coef(fitgf, terms = "time")
# Terms: cure probability model
coef(fitgf, terms = "cure")
# Object summaries
summary(fitgf)
# Information criterion
AIC(fitgf) # Akaike information criterion
BIC(fitgf) # Bayesian information criterion
# Extract Log-Likelihood
logLik(fitgf)
# Calculate Variance-Covariance Matrix for a Fitted Model Object
vcov(fitgf)

## fit Marshall-Olkin extended extreme value standart mixture model
fitmo <- curereg(Surv(FAILTIME, FAILCENS) ~ TRT + SEX + AGE, cureformula = ~ TRT + SEX + AGE,
                 data = e1684, timedist = "moeev", ncausedist = "bernoulli")

# Output of 'curereg' object
fitmo
# Extract Model Coefficients
coef(fitmo)
# Extract model coefficients:
# Terms: failure time distribution model
coef(fitmo, terms = "time")
# Terms: cure probability model
coef(fitmo, terms = "cure")
# Object summaries
summary(fitmo)
# Information criterion
AIC(fitmo) # Akaike information criterion
BIC(fitmo) # Bayesian information criterion
# Extract Log-Likelihood
logLik(fitmo)
# Calculate Variance-Covariance Matrix for a Fitted Model Object
vcov(fitmo)
}

\keyword{ models }
\keyword{ survival }
