-- Create tables mirroring the processed artifacts.
-- Suitable for SQLite; replace TEXT/INTEGER/VARCHAR sizes for Postgres if desired.

CREATE TABLE IF NOT EXISTS pubmed_gpt_staging (
  pmid TEXT NOT NULL,
  confidence TEXT NOT NULL CHECK (confidence IN ('YES','probablyYES','NO')),
  compound TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS new_inhibitors_processed (
  pmid TEXT NOT NULL,
  confidence TEXT NOT NULL CHECK (confidence IN ('YES','probablyYES')),
  compound TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS all_inhibitors_processed (
  compound TEXT NOT NULL,
  pubmed_references INTEGER NOT NULL CHECK (pubmed_references >= 0),
  known_status TEXT NOT NULL CHECK (known_status IN ('known','new'))
);

CREATE TABLE IF NOT EXISTS mesh_bioactive_staging (
  Type TEXT NOT NULL CHECK (Type IN ('Descriptor','SCR')),
  MeSH_UI TEXT NOT NULL,
  Name TEXT NOT NULL,
  OneTreeNumber TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS pubtator_filtered_staging (
  pmid TEXT NOT NULL,
  mention TEXT NOT NULL,
  normalized_id TEXT NOT NULL
);
