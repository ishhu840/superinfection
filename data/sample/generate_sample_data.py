"""Generate sample datasets for viral interference research."""

import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

def generate_sample_metadata(n_samples: int = 200) -> pd.DataFrame:
    """Generate sample metadata for mosquito samples.
    
    Args:
        n_samples: Number of samples to generate
        
    Returns:
        DataFrame with sample metadata
    """
    # Sample locations (major cities with Aedes aegypti presence)
    locations = [
        ('Miami', 'Florida', 'USA', 25.7617, -80.1918),
        ('Houston', 'Texas', 'USA', 29.7604, -95.3698),
        ('Bangkok', 'Bangkok', 'Thailand', 13.7563, 100.5018),
        ('Singapore', 'Singapore', 'Singapore', 1.3521, 103.8198),
        ('Rio de Janeiro', 'Rio de Janeiro', 'Brazil', -22.9068, -43.1729),
        ('Cairns', 'Queensland', 'Australia', -16.9186, 145.7781),
        ('Medellin', 'Antioquia', 'Colombia', 6.2442, -75.5812),
        ('Hanoi', 'Hanoi', 'Vietnam', 21.0285, 105.8542)
    ]
    
    # Generate sample data
    samples = []
    base_date = datetime(2020, 1, 1)
    
    for i in range(n_samples):
        # Random location
        city, state, country, lat, lon = locations[np.random.choice(len(locations))]
        
        # Random collection date (2020-2023)
        days_offset = np.random.randint(0, 1095)  # 3 years
        collection_date = base_date + timedelta(days=days_offset)
        
        # Season based on collection date
        month = collection_date.month
        if month in [12, 1, 2]:
            season = 'Winter'
        elif month in [3, 4, 5]:
            season = 'Spring'
        elif month in [6, 7, 8]:
            season = 'Summer'
        else:
            season = 'Fall'
        
        # Mosquito species (mostly Aedes aegypti)
        species = np.random.choice(
            ['Aedes_aegypti', 'Aedes_albopictus', 'Culex_quinquefasciatus'],
            p=[0.7, 0.2, 0.1]
        )
        
        # Collection method
        collection_method = np.random.choice(
            ['BG_trap', 'CDC_trap', 'Gravid_trap', 'Manual_collection'],
            p=[0.4, 0.3, 0.2, 0.1]
        )
        
        # Environmental factors
        temperature = np.random.normal(28, 5)  # Celsius
        humidity = np.random.normal(75, 15)  # Percentage
        rainfall = np.random.exponential(5)  # mm in past week
        
        # Pool information
        pool_size = np.random.randint(1, 51)  # 1-50 mosquitoes per pool
        
        samples.append({
            'sample_id': f'MOS_{i+1:04d}',
            'collection_date': collection_date.strftime('%Y-%m-%d'),
            'city': city,
            'state_province': state,
            'country': country,
            'latitude': lat + np.random.normal(0, 0.1),  # Add some noise
            'longitude': lon + np.random.normal(0, 0.1),
            'season': season,
            'mosquito_species': species,
            'collection_method': collection_method,
            'pool_size': pool_size,
            'temperature_c': round(temperature, 1),
            'humidity_percent': round(max(0, min(100, humidity)), 1),
            'rainfall_mm_week': round(max(0, rainfall), 1),
            'sequencing_depth': np.random.randint(1000000, 10000000),  # reads
            'host_reads_removed': np.random.randint(500000, 8000000)
        })
    
    return pd.DataFrame(samples)

def generate_viral_abundance_matrix(metadata: pd.DataFrame) -> pd.DataFrame:
    """Generate viral abundance matrix with realistic patterns.
    
    Args:
        metadata: Sample metadata
        
    Returns:
        DataFrame with viral abundances
    """
    n_samples = len(metadata)
    
    # Define virus groups with different prevalence patterns
    viruses = {
        # High prevalence viruses (common)
        'Aedes_aegypti_densovirus': {'base_prev': 0.6, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        'Culex_Y_virus': {'base_prev': 0.4, 'seasonal': False, 'species_specific': 'Culex_quinquefasciatus'},
        
        # Medium prevalence viruses
        'Phasi_Charoen_like_virus': {'base_prev': 0.3, 'seasonal': True, 'species_specific': None},
        'Hubei_mosquito_virus_1': {'base_prev': 0.25, 'seasonal': False, 'species_specific': None},
        'Aedes_flavivirus': {'base_prev': 0.2, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        
        # Lower prevalence viruses
        'Menghai_rhabdovirus': {'base_prev': 0.15, 'seasonal': False, 'species_specific': None},
        'Xinzhou_spider_virus': {'base_prev': 0.12, 'seasonal': True, 'species_specific': None},
        'Aedes_anphevirus': {'base_prev': 0.1, 'seasonal': False, 'species_specific': 'Aedes_aegypti'},
        
        # Rare viruses
        'Chikungunya_virus': {'base_prev': 0.05, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        'Dengue_virus': {'base_prev': 0.03, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        'Zika_virus': {'base_prev': 0.02, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        'Yellow_fever_virus': {'base_prev': 0.01, 'seasonal': True, 'species_specific': 'Aedes_aegypti'},
        
        # Additional insect-specific viruses
        'Lammi_virus': {'base_prev': 0.08, 'seasonal': False, 'species_specific': None},
        'Aedes_alboannulatus_virus': {'base_prev': 0.06, 'seasonal': True, 'species_specific': 'Aedes_albopictus'},
        'Culex_tritaeniorhynchus_rhabdovirus': {'base_prev': 0.04, 'seasonal': False, 'species_specific': 'Culex_quinquefasciatus'}
    }
    
    # Initialize abundance matrix
    abundance_matrix = pd.DataFrame(
        0.0, 
        index=metadata['sample_id'], 
        columns=list(viruses.keys())
    )
    
    # Generate abundances for each virus
    for virus, params in viruses.items():
        base_prevalence = params['base_prev']
        
        for idx, row in metadata.iterrows():
            sample_id = row['sample_id']
            
            # Adjust prevalence based on factors
            prevalence = base_prevalence
            
            # Species specificity
            if params['species_specific'] and row['mosquito_species'] != params['species_specific']:
                prevalence *= 0.1  # Much lower in other species
            
            # Seasonal effects
            if params['seasonal']:
                if row['season'] in ['Summer', 'Fall']:
                    prevalence *= 1.5  # Higher in warmer seasons
                else:
                    prevalence *= 0.7
            
            # Temperature effects
            temp_effect = 1 + (row['temperature_c'] - 25) * 0.02
            prevalence *= max(0.1, temp_effect)
            
            # Humidity effects (some viruses prefer higher humidity)
            if virus in ['Dengue_virus', 'Zika_virus', 'Chikungunya_virus']:
                humidity_effect = 1 + (row['humidity_percent'] - 70) * 0.01
                prevalence *= max(0.5, humidity_effect)
            
            # Pool size effect (larger pools more likely to have virus)
            pool_effect = 1 + (row['pool_size'] - 25) * 0.005
            prevalence *= pool_effect
            
            # Determine if virus is present
            if np.random.random() < min(0.95, prevalence):
                # Generate abundance (log-normal distribution)
                abundance = np.random.lognormal(mean=8, sigma=2)
                abundance_matrix.loc[sample_id, virus] = abundance
    
    # Add some viral exclusion patterns
    # Dengue and Zika rarely co-occur
    dengue_positive = abundance_matrix['Dengue_virus'] > 0
    zika_positive = abundance_matrix['Zika_virus'] > 0
    coinfected = dengue_positive & zika_positive
    
    # Remove 80% of co-infections
    remove_coinfection = np.random.random(sum(coinfected)) < 0.8
    coinfected_samples = abundance_matrix.index[coinfected]
    samples_to_clear = coinfected_samples[remove_coinfection]
    
    for sample in samples_to_clear:
        # Randomly choose which virus to remove
        if np.random.random() < 0.5:
            abundance_matrix.loc[sample, 'Dengue_virus'] = 0
        else:
            abundance_matrix.loc[sample, 'Zika_virus'] = 0
    
    # Chikungunya and Dengue also show some exclusion
    chik_positive = abundance_matrix['Chikungunya_virus'] > 0
    dengue_positive = abundance_matrix['Dengue_virus'] > 0
    coinfected = chik_positive & dengue_positive
    
    remove_coinfection = np.random.random(sum(coinfected)) < 0.6
    coinfected_samples = abundance_matrix.index[coinfected]
    samples_to_clear = coinfected_samples[remove_coinfection]
    
    for sample in samples_to_clear:
        if np.random.random() < 0.5:
            abundance_matrix.loc[sample, 'Chikungunya_virus'] = 0
        else:
            abundance_matrix.loc[sample, 'Dengue_virus'] = 0
    
    # Some viruses show positive associations
    # Aedes-specific viruses tend to co-occur
    aedes_viruses = ['Aedes_aegypti_densovirus', 'Aedes_flavivirus', 'Aedes_anphevirus']
    
    for i, row in metadata.iterrows():
        sample_id = row['sample_id']
        if row['mosquito_species'] == 'Aedes_aegypti':
            # If one Aedes virus is present, increase chance of others
            aedes_present = [v for v in aedes_viruses if abundance_matrix.loc[sample_id, v] > 0]
            
            if len(aedes_present) > 0:
                for virus in aedes_viruses:
                    if abundance_matrix.loc[sample_id, virus] == 0 and np.random.random() < 0.3:
                        abundance = np.random.lognormal(mean=7, sigma=1.5)
                        abundance_matrix.loc[sample_id, virus] = abundance
    
    return abundance_matrix

def generate_sample_datasets():
    """Generate and save sample datasets."""
    print("Generating sample datasets for viral interference research...")
    
    # Create output directory
    output_dir = Path('/Users/ishtiaq/Desktop/super-virus/data/sample')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate metadata
    print("Generating sample metadata...")
    metadata = generate_sample_metadata(n_samples=200)
    metadata.to_csv(output_dir / 'sample_metadata.csv', index=False)
    print(f"Saved metadata for {len(metadata)} samples")
    
    # Generate viral abundance matrix
    print("Generating viral abundance matrix...")
    abundance_matrix = generate_viral_abundance_matrix(metadata)
    abundance_matrix.to_csv(output_dir / 'viral_abundance_matrix.csv')
    print(f"Saved abundance matrix: {abundance_matrix.shape[0]} samples x {abundance_matrix.shape[1]} viruses")
    
    # Generate summary statistics
    print("\nDataset Summary:")
    print(f"Total samples: {len(metadata)}")
    print(f"Date range: {metadata['collection_date'].min()} to {metadata['collection_date'].max()}")
    print(f"Countries: {', '.join(metadata['country'].unique())}")
    print(f"Mosquito species: {', '.join(metadata['mosquito_species'].unique())}")
    
    print(f"\nViral prevalence:")
    prevalence = (abundance_matrix > 0).mean().sort_values(ascending=False)
    for virus, prev in prevalence.items():
        print(f"  {virus}: {prev:.1%}")
    
    print(f"\nSample viral richness:")
    richness = (abundance_matrix > 0).sum(axis=1)
    print(f"  Mean: {richness.mean():.1f} viruses per sample")
    print(f"  Range: {richness.min()}-{richness.max()} viruses per sample")
    
    # Save summary
    summary = {
        'total_samples': len(metadata),
        'total_viruses': len(abundance_matrix.columns),
        'date_range': f"{metadata['collection_date'].min()} to {metadata['collection_date'].max()}",
        'countries': metadata['country'].unique().tolist(),
        'species': metadata['mosquito_species'].unique().tolist(),
        'viral_prevalence': prevalence.to_dict(),
        'mean_viral_richness': float(richness.mean()),
        'viral_richness_range': [int(richness.min()), int(richness.max())]
    }
    
    import json
    with open(output_dir / 'dataset_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nSample datasets generated successfully!")
    print(f"Files saved to: {output_dir}")
    print("- sample_metadata.csv")
    print("- viral_abundance_matrix.csv")
    print("- dataset_summary.json")

if __name__ == '__main__':
    generate_sample_datasets()