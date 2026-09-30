# Mejoras pendientes – PCR_Virtual.py

Lista surgida de la revisión del script (2026-09-30). Los puntos 1 y 3-5 fueron confirmados con pruebas.

Prioridad sugerida: 1 y 2 primero (afectan qué se encuentra); luego 3-5 y 7-8 (arreglos rápidos); luego 10, 11 y 13 (valor biológico).

## 🔴 Fallos importantes

- [ ] **1. Solo se detecta el mejor sitio de unión por secuencia.** `aligner.align()` devuelve solo los alineamientos de puntaje máximo, y el filtro `-p` se aplica sobre ellos. Sitios con menor puntaje que superan `-p` se pierden. Ej.: sitio perfecto + sitio con 1 mismatch (0.95) en el mismo registro → solo sale el perfecto. Si el partidor une perfecto en una zona irrelevante y con 1 mismatch en el blanco real, el producto real no aparece.
- [ ] **2. Bases degeneradas (IUPAC) cuentan como mismatch.** Con `CGTCCAARRGGATACTGATC` el máximo posible es 18/20 = 0.90 (por eso el ejemplo usa `-p 0.85`). Usar una matriz de sustitución que acepte R, Y, N, etc. Lo mismo para `N` en el genoma.
- [ ] **3. `-s` con valor inválido genera FASTA roto** (encabezados sin secuencia). Usar `choices=['A', 'a']`.
- [ ] **4. `-l` mal escrito provoca error** (`-l 100` → `IndexError`). Validar formato y que min ≤ max.
- [ ] **5. `-p` como porcentaje devuelve vacío sin avisar** (`-p 85` no encuentra nada; espera 0.85). El help dice "Porcentaje". Validar rango 0-1 o aceptar ambos formatos.

## 🟡 Fallos menores / robustez

- [ ] **6. Coordenadas negativas posibles** en uniones parciales al inicio de la secuencia (`c_1 = subali_t[0] - dif_L`). Un slice negativo toma desde el final → secuencia del partidor vacía o incorrecta.
- [ ] **7. Códigos de salida incorrectos.** `exit()` al faltar Biopython o argumentos devuelve 0. Usar `sys.exit(1)`, o `required=True` / `parser.error()`.
- [ ] **8. `-r` y `-rc` juntos:** `-rc` sobrescribe a `-r` sin avisar. Hacerlos mutuamente excluyentes.
- [ ] **9. Partidores sin validar** (caracteres no nucleotídicos no dan error).

## 🔵 Mejoras biológicas

- [ ] **10. Mismatches en el extremo 3′** deberían pesar más (un mismatch en la última base suele impedir la amplificación). Opción tipo "máximo N mismatches en las últimas X bases".
- [ ] **11. Genomas circulares** (ej. HPV): no se detectan amplicones que crucen la posición 1. Flag `--circular` que concatene el inicio de la secuencia al final.
- [ ] **12. Productos anidados o redundantes:** se combinan todos los forward con todos los reverse. Máximo de `-l` por defecto más realista (ej. 5 kb) o reportar solo el producto más corto por forward.

## 🟢 Salida y usabilidad

- [ ] **13. Coordenadas y hebra en el encabezado**, ej. `>ID:1234-1653(+) 409pb`. En los productos `RC_`, convertir las coordenadas a la hebra original.
- [ ] **14. Resumen final:** `resume` se llena pero no se usa. Imprimir en stderr productos por secuencia y secuencias sin producto.
- [ ] **15. Opción `-o`** para escribir a un archivo.
- [ ] **16. Help:** dice `PCR_VIRTUAL.py` (el archivo es `PCR_Virtual.py`). La lógica `-r`/`-rc` es poco intuitiva: el reverse tal como se encarga (5′→3′) va en `-rc`. Aclarar o invertir nombres.
- [ ] **17. Mayúsculas/minúsculas con `-c`:** la convención está invertida respecto al modo normal (partidores en mayúscula). Evaluar unificar.

## ⚪ Calidad de código

- [ ] **18. Código duplicado** entre hebra + y − (~30 líneas). Extraer a una función `buscar_productos(seq, prefijo)`.
- [ ] **19. Cálculos repetidos:** `reverse_complement()` se calcula 3-4 veces por registro y el `PairwiseAligner` se crea en cada llamada.
- [ ] **20. Código sin efecto:** `type(forward and reverse) == list` siempre es verdadero; bloque comentado muerto; condiciones `> c_1` / `> c_2` sobran.
- [ ] **21. `args` global:** mover el parseo a `main()` y pasar parámetros explícitos para poder importar y testear. El docstring de `PCR()` dice `str`, pero recibe `Seq`.
- [ ] **22. Sin tests y README mínimo.** Agregar casos de prueba y ejemplos de uso.
