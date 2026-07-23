# BibliotecaMatrices

Biblioteca matricial inmutable para Java 21 con elementos
`ch.obermuhlner.math.big.BigComplex`, precisión decimal configurable y dos
representaciones independientes:

- `array.Matriz` almacena realmente `BigComplex[][]`.
- `lista.Matriz` almacena realmente `ArrayList<ArrayList<BigComplex>>`.

`lista.Matriz` no construye, llama ni delega en `array.Matriz`. Los algoritmos
compartidos de `comun.AbstractMatriz` solo usan la interfaz de acceso
`MatrizCompleja`; el método factoría abstracto `crear(...)` garantiza que cada
resultado conserva la representación de su receptor.

## Requisitos y ejecución

- JDK 21.
- Maven 3.9 o posterior.
- Docker, solo para la ejecución contenerizada.

```bash
mvn clean test
mvn exec:java
docker build -t biblioteca-matrices .
docker run --rm biblioteca-matrices
```

La prueba `CoberturaApiTest` ejecuta `BibliotecaMatricesApp`, descubre por reflexión
los constructores y métodos públicos, los invoca individualmente y falla ante
cualquier método no probado, retorno nulo o excepción inesperada. El principal
imprime entradas, resultados, `CORRECTO` y un resumen de cobertura.

## Dependencias

- `ch.obermuhlner:big-math:2.3.2`, que proporciona `BigComplex`.
- `org.junit.jupiter:junit-jupiter:5.11.4`, solo para pruebas.

No se almacena ni se convierte una matriz a `double` para ejecutar algoritmos.
Los constructores heredados con `double` convierten usando la representación
decimal de `Double.toString()`.

## Estructura

```text
src/main/java/
├── BibliotecaMatricesApp.java
├── array/Matriz.java
├── lista/Matriz.java
├── comun/
│   ├── AbstractMatriz.java
│   ├── Complejos.java
│   ├── MatrizCompleja.java
│   ├── excepciones/
│   └── resultados/
└── otros/
    ├── MatrizBooleana.java
    └── MatrizModular.java
```

`otros` separa correctamente el semianillo booleano y los cuerpos finitos
primos de la aritmética compleja.

## Precisión, tolerancia y mutabilidad

La política predeterminada es `MathContext.DECIMAL128` y tolerancia `1E-28`.
Cada matriz conserva ambas propiedades de forma inmutable. Una operación entre
matrices con políticas distintas se rechaza para no mezclar precisiones
silenciosamente. Los algoritmos redondean explícitamente con el contexto del
receptor.

Se rechazan matrices con cero filas o cero columnas. Todos los constructores,
getters, conversiones y resultados copian los contenedores; ninguna estructura
interna queda expuesta. `BigComplex` es inmutable.

La igualdad `equals` es exacta según la igualdad matemática de `BigComplex`;
`igualExacta` la hace explícita y `aproximadamenteIgual` recibe tolerancia.

## Catálogo de API

Las dos implementaciones ofrecen la misma API de instancia:

- Construcción y acceso: dimensiones, filas, columnas, elementos, copias,
  conversión entre representaciones, submatrices, inserción, eliminación,
  concatenación, redimensión, `map` y `zip`.
- Aritmética: suma, resta, opuesta, escalares, productos matricial, Hadamard,
  Kronecker, Khatri-Rao, interior, exterior y Frobenius, potencias, formas
  bilineal, sesquilineal y cuadrática.
- Simetrías: transpuesta, conjugada, adjunta hermítica, partes real,
  imaginaria, simétrica, antisimétrica, hermítica y antihermítica, diagonales y
  partes triangulares.
- Clasificación: estructura diagonal, triangular, banda, Hessenberg,
  Toeplitz, Hankel, circulante, permutación, monomial y por bloques;
  propiedades estocásticas, simétricas, hermíticas, normales, unitarias,
  nilpotentes, periódicas, de rango, definitud, proyección, rotación,
  reflexión, reducibilidad y dominancia diagonal.
- Eliminación: operaciones elementales de filas y columnas, Gauss,
  Gauss-Jordan, RREF, pivotes, rango y nulidad, con traza detallada.
- Determinantes: Gauss, Laplace, LU, Bareiss compatible, 2×2, Sarrus,
  menores, cofactores, adjunta clásica y Schur.
- Inversas y sistemas: Gauss-Jordan, adjunta, LU, diagonal, triangular,
  bloques, pseudoinversa, ecuaciones matriciales, Rouché-Frobenius, Cramer,
  mínimos cuadrados y bases de los cuatro espacios fundamentales.
- Normas y condición: traza, normas uno, infinito, máximo, Frobenius,
  espectral, nuclear y dos; radio espectral, condición y ángulo Frobenius.
- Espectro: polinomio característico por Faddeev-LeVerrier, autovalores por
  iteración QR, autopares y residuos, espacios propios, cadenas de Jordan,
  Cayley-Hamilton, potencia y cociente de Rayleigh.
- Formas y factorizaciones: diagonalización general, ortogonal y unitaria,
  inercia, LU/LUP/LDU, QR clásico/modificado/Householder/Givens, Cholesky,
  LDL*, SVD, polar, Schur, Hessenberg, bidiagonal, rango, CUR y NMF.
- Ortogonalidad: normalización, Gram-Schmidt, Gram, bases ortonormales,
  complementos y proyectores.
- Funciones: polinomios, exponencial por serie/escalado y diagonalización,
  logaritmo y raíz principales, seno, coseno, tangente, hiperbólicas, signo,
  valor absoluto y funciones racionales o analíticas.
- Cálculo: diferenciales de producto, inversa y determinante, gradientes,
  jacobiano, hessiano, derivación central e integración trapezoidal.
- Iterativos: Jacobi, Gauss-Seidel, SOR, gradiente conjugado, GMRES, BiCG,
  BiCGSTAB, potencia, potencia inversa, Rayleigh, QR, Lanczos, Arnoldi,
  refinamiento y precondicionadores.
- Bloques, estadística y aplicaciones: partición, ensamblado, Schur,
  covarianza, correlación, PCA, blanqueamiento, regresión, matriz sombrero,
  confusión y compresión SVD.
- Grafos y Markov: grados, laplacianas, recorridos, clausura transitiva,
  componentes, conectividad, distribuciones tras pasos, estacionaria y
  regularidad.
- Vectorización: orden por filas o columnas, `reshape`, `flatten`,
  desvectorización, suma de Kronecker y la identidad `vec(AXB)`.
- Verificación: inversa, LU, QR, Cholesky, SVD, diagonalización, Schur, Jordan,
  pseudoinversa, sistema, autopar, proyector, unitaria, hermítica,
  Cayley-Hamilton y residuos.

Las factorías de ambas clases incluyen cero, identidades rectangulares,
diagonal, escalar, fila, columna, triangulares, bi/tridiagonales, banda,
Hessenberg, Toeplitz, Hankel, circulante, Vandermonde, Hilbert, Fourier,
Hadamard, permutación, rotaciones 2D/3D, reflexión, proyección, Householder,
Givens, adyacencia, incidencia, compañera, bloque de Jordan y matrices
aleatorias normales, dispersas, hermíticas, unitarias y definidas positivas.

## Convenciones matemáticas

- Los índices empiezan en cero; los intervalos de submatriz son
  `[desde, hasta)`.
- El producto interior es sesquilineal y conjuga el primer argumento.
- `adjuntaClasica` es la transpuesta de cofactores;
  `adjuntaHermitica` es la transpuesta conjugada.
- La definitud compleja exige hermiticidad y se interpreta mediante `x*Ax`.
- La pseudoinversa verifica las cuatro ecuaciones de Moore-Penrose.
- La SVD es reducida; `svdTruncada(r)` conserva los `r` mayores valores.
- La diagonalización ortogonal exige simetría real y la unitaria normalidad.
- Las funciones «principal» usan la rama principal y rechazan hipótesis que no
  permiten certificarla con la precisión vigente.
- Una matriz de transición es cuadrada, no negativa y estocástica por filas.
- La conectividad de un grafo dirigido se evalúa sobre su grafo subyacente no
  dirigido.

## Complejidad aproximada

| Operación | Tiempo | Memoria adicional |
|---|---:|---:|
| suma / Hadamard | `O(mn)` | `O(mn)` |
| producto | `O(mnp)` | `O(mp)` |
| Gauss, LU, inversa | `O(n³)` | `O(n²)` |
| QR | `O(mn²)` | `O(mn)` |
| Cholesky | `O(n³)` | `O(n²)` |
| iteración QR espectral | `O(k n³)` | `O(n²)` |
| SVD espectral | `O(n³ + mn²)` | `O(mn)` |
| Kronecker | `O(mnpq)` | tamaño del resultado |
| clausura transitiva | `O(n³)` | `O(n²)` |

## Limitaciones matemáticas genuinas

Los autovalores generales son numéricos, no expresiones radicales simbólicas.
El resultado depende de la tolerancia y del máximo de iteraciones. Una matriz
puede no ser diagonalizable; una raíz o logaritmo matricial puede no existir o
no ser único; Cholesky exige hermiticidad definida positiva; Cramer exige
cuadrada regular; la matriz de Hadamard de Sylvester exige orden potencia de
dos; la inversa modular solo se ofrece sobre módulos primos. Estas situaciones
se comunican con excepciones de dominio y no con valores ficticios.
