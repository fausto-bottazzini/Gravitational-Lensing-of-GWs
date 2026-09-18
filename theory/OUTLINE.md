# `theory/` — guía de la carpeta

## Qué es `theory.tex`

Un texto autocontenido, escrito como libro de texto breve y no como informe,
sobre lente gravitacional de ondas gravitacionales en régimen de óptica de
ondas, especializado en el caso en que la fuente y el lente pertenecen al mismo
triple jerárquico.

La dependencia con el resto del proyecto va en un solo sentido: `report/` cita
a este texto, y este texto se lee sin `report/`. No contiene rutas del
repositorio, nombres de tests ni narrativa de implementación; los valores
numéricos concretos aparecen únicamente en recuadros `ejemplo`, separados del
hilo teórico.

## Estructura

Siete capítulos, un apéndice de deducciones y un prefacio. 41 páginas, con
una única figura (4.1, el tren de pulsos, en §4.3).

| | Capítulo | Contenido |
|---|---|---|
| 1 | Triples jerárquicos como sistemas de lente | Anatomía y notación del sistema; la condición de jerarquía desglosada en sus cuatro exigencias físicas (estabilidad de Mardling–Aarseth, radio de Hill, Kozai–Lidov y su supresión relativista, estacionariedad de la órbita externa); la binaria interna como fuente puntual; qué cambia en las distancias cuando el lente no es cosmológico; los dos regímenes temporales; las escalas del problema |
| 2 | De la ecuación de onda a la integral de difracción | Propagación de una GW sobre fondo curvo; reducción a campo escalar en el límite de longitud de onda corta; ecuación de Helmholtz e índice de refracción; lente delgado e integral de Fresnel–Kirchhoff; forma adimensional en $(w,y)$, convención de Fourier y fase de referencia |
| 3 | El lente de masa puntual | Potencial $\psi=\ln\lvert x\rvert$; óptica geométrica (ecuación de lente, magnificaciones, Paczyński); retardo entre imágenes y $F_\text{geo}$ con índice de Morse; solución exacta hipergeométrica y sus tres límites; forma radial equivalente con Bessel |
| 4 | Geometría del triple como sistema de lente | Órbita externa proyectada; $D_{LS}$ como proyección sobre la línea de visión y el radio de Einstein dependiente del tiempo; lente repetido; cinemática (Roemer, Doppler, corrimiento orbital); límites de la aproximación paraxial y precisión global del tratamiento |
| 5 | Más allá de la aproximación de campo débil | Propagación tensorial (Regge–Wheeler, Zerilli, Teukolsky); fenómenos que el tratamiento escalar no puede producir; por qué queda fuera de alcance |
| 6 | La binaria interna como fuente | Inspiral cuasi-circular y masa de chirp; forma de onda en frecuencia por fase estacionaria (TaylorF2 a 2PN); modelos IMR |
| 7 | Los dos regímenes | Lente estático (multiplicación en frecuencia) y lente en movimiento (aproximación adiabática); cuadratura entre el pulso de lente y la modulación Doppler |
| A | Deducciones | Diez deducciones paso a paso: ecuación de propagación sobre fondo curvo; Helmholtz; Fresnel–Kirchhoff; forma adimensional; $\psi$ de una masa puntual; jacobiano axisimétrico; fase estacionaria; reducción radial; la forma cerrada del lente puntual; masa de chirp |

Convención de enlaces: cada resultado enunciado en el cuerpo remite a su
deducción completa en el apéndice, y cada sección del apéndice está referenciada
desde la ecuación correspondiente. Ese es el mecanismo que permite que el cuerpo
se lea corrido sin perder las cuentas.

## Convenciones fijadas en el texto

- Signatura $(-,+,+,+)$.
- $G$ y $c$ explícitas en todas las fórmulas, sin excepción, para que cada
  expresión sea evaluable en SI tal como está escrita.
- $f$ es siempre la frecuencia de la onda gravitacional, el doble de la orbital.
- Transformada de Fourier con $h(f)=\int h(t)\,e^{-2\pi ift}\,dt$. La literatura
  de lentes suele usar el signo opuesto, y por eso las expresiones de $F$
  aparecen conjugadas respecto de las de Takahashi & Nakamura (2003).
- Notación: $h(f)$ sin lente, $h_L(f)=F(f)h(f)$ a la salida del lente,
  $F(w,y)$ factor de amplificación, $w=8\pi GM_Lf/c^3$ frecuencia adimensional,
  $y$ parámetro de impacto en unidades del radio de Einstein.

## Bibliografía

`refs.bib` es compartido con `report/report.tex` (que lo incluye como
`../theory/refs`). Los PDF no están versionados, porque son artículos con
copyright de sus revistas: `bibliography/README.md` lista el DOI de cada uno
y el enlace a arXiv de los ocho que tienen preprint libre, que es lo que hace
falta para descargarlos. Las citas indican la ecuación o sección exacta del
trabajo citado cuando se toma de él un resultado concreto; las cuatro
referencias que son libros (Maggiore; Schneider, Ehlers & Falco; Murray &
Dermott; Chandrasekhar) se citan por capítulo y sección.

## Archivos

| Archivo | Para qué |
|---|---|
| `theory.tex` | fuente del documento, **y la manera de leerlo si sos un agente** (ver abajo) |
| `theory.pdf` | documento compilado, que es lo que se entrega |
| `refs.bib` | bibliografía, compartida con `report/` |
| `pulsos_esquema.py` | genera la figura 4.1 (el tren de pulsos, §4.3), esquemática y sin números: `theory.pdf` enuncia el fenómeno, no la instancia — los números de este proyecto están en `cases/case_B_monochromatic/caseB_repeated_pulses.png` |
| `pulsos_esquema.pdf` | salida del script anterior, incluida por `theory.tex` |
| `paraxial_validity.py` | genera los números del Cuadro 4.1 (validez paraxial) a partir de `src/gwlens/system.py`; lo corre `reproduce.sh` |
| `paraxial_validity_numbers.json` | salida del script anterior, que es lo que el Cuadro 4.1 transcribe |
| `OUTLINE.md` | este archivo |

Los archivos auxiliares de LaTeX (`.aux`, `.log`, `.out`, `.toc`, `.bbl`,
`.blg`) no forman parte del repositorio y se regeneran al compilar.

## Para leerlo desde un agente: `theory.tex`, no `theory.pdf`

Leer el PDF cuesta caro, porque se lee renderizando cada página como imagen.
**La fuente es la lectura barata y es la misma información**: unos 135 000
caracteres de `.tex` contra un PDF de unos 540 KB en imágenes. El `.tex` se lee
corrido sin problema — los acentos van como `\'a`, la matemática entre `$`, y
la estructura (`\chapter`, `\section`, `\begin{ejemplo}`) es la del índice de
más arriba.

Tampoco hace falta generar una versión Markdown. Existió una
(`theory.md`, borrada el 2026-09-17) que se justificaba como "la versión
económica en tokens": medida en su momento, eran 126 915 caracteres
contra 133 091 del
`.tex`, o sea un 5 % menos, a cambio de mantenerla sincronizada a mano — y ya
se había desincronizado una vez, quedando un archivo que decía "generado desde
theory.tex" y no lo estaba, que es peor que no tenerlo. Un volcado línea a
línea no es un resumen: no aporta nada que no aporte la fuente. Si alguna vez
hace falta algo más barato que la fuente, tiene que ser un resumen de verdad
—resultados, condiciones de validez, dónde está cada deducción— y no una
conversión automática.

## Compilación

Tres pasadas de `pdflatex` con `bibtex` entre la primera y la segunda, que es lo
que hace `reproduce.sh`:

```sh
cd theory
pdflatex -interaction=nonstopmode theory.tex
bibtex theory
pdflatex -interaction=nonstopmode theory.tex
pdflatex -interaction=nonstopmode theory.tex
```

Requiere `newtx`, `babel-spanish`, `microtype`, `titlesec`, `fancyhdr`,
`tcolorbox`, `tikz`, `physics`, `siunitx` y `natbib`.
