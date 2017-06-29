# fit_walk

Cómo funciona la simulación.
La función principal recibe como argumentos:
$n=10$, 			Número de dimensiones
$tp=1$, 			tamaño del paso
$N=1000$, 		tamaño de la población
$max.time=10000$, 	número máximo de mutaciones (iteraciones)
$amb.tau=NA$, 		la periodicidad con que ocurre un cambio ambiental (se alimenta como 1/tau)
$var.amb$=1		Amplitud de la variabilidad ambiental
(se indican los valores que las variables asumen por defecto)

Algoritmo:

1) Se determina la posición (pos) inicial (n dimensiones) en 0.5
2) Se determina la posición del óptimo (optimo)
3) Antes de en entrar en el ciclo principal se determina el tiempo de ocurrencia (ciclo) del primer evento de variabilidad ambiental a partir de una distribución exponencial con promedio en amb.tau (de no consiserar un escenario sin variabilidad ambiental (amb.tau=NA))
4) Se inicia el ciclo con un while que concluye cuando los ciclos igualan a max.time.
5) Se controla si el ciclo en curso corresponde al momento de un evento de var.amb. en cuyo caso el óptimo se desplaza optimo <- o.muller.v(n,rnorm(1,mean=0,sd=var.amb))). Además se calcula el momento del sisuiente evento de variabilidad ambiental (como en 3).
6) Se calcula un vector (mut) isotrópico según algoritmo de muller  en n dimensiones que representa la mutación
7) Se calcula la nueva posición (new.pos), suma vectorial mut + pos.
8) Se calcula Zwt = sqrt(sum((pos-optimo)^2))
9) Se calcula Zmut <- sqrt(sum((new.pos-optimo)^2))
10) Una variable true/false (advantageos) marca true si la nueva posición está más cerca del óptimo que en el paso anterior.
11) Se calcula s = exp((Zwt^2-Zmut^2)/2)-1
12) Se calculan las probabilidades de fijación p <- p.fix(s,N) y p.Ninf <- p.fix(s,1e5). Según la fórmula p.fix <- function(s,N) (1-exp(-2*s))/(1-exp(-2*N*s))
13) Se calculan las components para el cálculo de F (suma de dos integrales) f.drift <- f.drift + (p-p.Ninf) y  f.sel <- f.sel + p.Ninf
14) Se determina si la nueva mutación se fija o no en un sorteo considerando p (alculado en 12). Si el sorteo sale positivo, se calcula la distancia x que se avanza con la nueva mutación fijada, la distancia se escala para comparar entre diferentes números de dimensiones por sqrt(n). La distancia se almacena en un vector de largo variable R.adv o R.del según si la fijación fue ventajosa o deletérea (valor de advantageous). Luego la variable pos, que contiene el valor actual de la posición en la caminata asume el valor de new.pos. El ciclo comienza devuevo.
15) La función devuelve una lista list(fixs=muts,adv=R.adv,del=R.del,drift=f.drift,selection=f.sel,N=N,muts=max.time,K=muts/max.time). Contiene (en orden) el número de fijaciones, los vectores de mutaciones  ventajosas y deletéreas, las componentes integradas de drift y selección, el tamaño poblacional, las mutaciones totales (ciclos de la simulación) y K que es fixs/muts.
Por ejemplo hoy chequeamos a partir de esta lista de multiples y versátiles datos, que K*N es igual a N*(drift+selection)/muts. Ya que ambos son F según las definiciones.
Otras funciones ad-hoc pueden llamar a esta función principal para barrer diferentes valores de uno o varios parámetros, colectar los resultados y construir la(s) curva(s) que se desee….

Chan….. chan…!!!!
Escribir código es entrete pero escribir su explicación es lata pura!!!!!! Jajaja.
