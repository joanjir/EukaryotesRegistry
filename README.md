# 🧬 EukaryotesRegistry  
### Sistema de consulta, evaluación y clasificación de genomas eucariontes desde la API de NCBI Datasets

---

## 📘 Descripción general

**EukaryotesRegistry** es un sistema automatizado para **consultar, filtrar y evaluar genomas eucariontes** disponibles en la base de datos **NCBI Datasets**, priorizando aquellos con **mayor calidad de ensamblaje**.

El objetivo principal es construir un registro de especies representativas por clase y filo, seleccionando los **mejores ensamblajes genómicos** disponibles (RefSeq o GenBank).  
El sistema realiza consultas jerárquicas (`Phylum → Class → Species`) y analiza métricas de calidad estructural de los genomas (cobertura, nivel de ensamblaje, N50 de contigs y scaffolds, número de scaffolds, etc.), con ejecución paralelizada y registro detallado de eventos.

---

## 🧩 Arquitectura del proyecto

EukaryotesRegistry/
│
├── config.py          # Configuración general, logging y control de concurrencia
├── http_client.py     # Cliente HTTP con reintentos y control de peticiones
├── taxonomy.py        # Resolución de taxID, clases y especies mediante NCBI Taxonomy API
├── genomes.py         # Evaluación y puntuación de ensamblajes genómicos
├── pipeline.py        # Flujo principal (consulta paralelizada y selección de mejores genomas)
├── main.py            # Punto de entrada del sistema
├── README.md          # Documentación general del proyecto
└── pipeline.log       # Registro de ejecución detallado
