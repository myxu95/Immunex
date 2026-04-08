import tempfile
import unittest
from pathlib import Path

import pandas as pd

from immunex.core.context import PipelineContext
from immunex.pipeline.nodes.contact_heatmap_node import ContactHeatmapNode


class TestContactHeatmapNode(unittest.TestCase):
    def test_generate_global_and_region_heatmaps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_path = Path(tmpdir)
            contacts_dir = temp_path / 'results' / 'analysis' / 'contacts'
            cdr1_dir = contacts_dir / 'cdr1'
            cdr1_tables = cdr1_dir / 'tables'
            cdr1_heatmaps = cdr1_dir / 'heatmaps'
            cdr1_tables.mkdir(parents=True)
            cdr1_heatmaps.mkdir(parents=True)

            report_df = pd.DataFrame(
                [
                    {
                        'interaction_class': 'peptide_tcr',
                        'phla_region': 'peptide',
                        'phla_chain_id': 'C',
                        'phla_resid': 5,
                        'phla_resname': 'TYR',
                        'phla_residue': 'TYR5',
                        'tcr_chain_id': 'D',
                        'tcr_chain': 'alpha',
                        'tcr_resid': 31,
                        'tcr_resname': 'SER',
                        'tcr_residue': 'SER31',
                        'tcr_region': 'CDR1',
                        'tcr_region_detailed': 'CDR1_alpha',
                        'contact_frequency': 0.8,
                        'contact_frames': 8,
                        'total_frames': 10,
                        'min_distance_angstrom': 3.2,
                    },
                    {
                        'interaction_class': 'hla_tcr',
                        'phla_region': 'HLA_alpha',
                        'mhc_region': 'groove',
                        'phla_chain_id': 'A',
                        'phla_resid': 68,
                        'phla_resname': 'LYS',
                        'phla_residue': 'LYS68',
                        'tcr_chain_id': 'E',
                        'tcr_chain': 'beta',
                        'tcr_resid': 95,
                        'tcr_resname': 'TRP',
                        'tcr_residue': 'TRP95',
                        'tcr_region': 'CDR3',
                        'tcr_region_detailed': 'CDR3_beta',
                        'contact_frequency': 0.6,
                        'contact_frames': 6,
                        'total_frames': 10,
                        'min_distance_angstrom': 3.5,
                    },
                ]
            )
            report_file = contacts_dir / 'contact_report.csv'
            contacts_dir.mkdir(parents=True, exist_ok=True)
            report_df.to_csv(report_file, index=False)

            cdr1_df = report_df[report_df['tcr_region'] == 'CDR1'].copy()
            cdr1_contacts_file = cdr1_tables / 'contacts.csv'
            cdr1_df.to_csv(cdr1_contacts_file, index=False)

            cdr2_contacts_file = cdr1_tables / 'empty_cdr2.csv'
            cdr3_contacts_file = cdr1_tables / 'empty_cdr3.csv'
            non_cdr_contacts_file = cdr1_tables / 'empty_non_cdr.csv'
            empty_df = pd.DataFrame(columns=report_df.columns)
            empty_df.to_csv(cdr2_contacts_file, index=False)
            empty_df.to_csv(cdr3_contacts_file, index=False)
            empty_df.to_csv(non_cdr_contacts_file, index=False)

            context = PipelineContext(
                system_id='demo',
                topology=str(temp_path / 'md.tpr'),
                trajectory_raw=str(temp_path / 'md.xtc'),
                output_dir=str(temp_path / 'results'),
            )
            context.results['contact_annotation'] = {
                'contact_report_file': str(report_file),
                'cdr1': {
                    'contacts_file': str(cdr1_contacts_file),
                    'heatmaps_dir': str(cdr1_heatmaps),
                },
                'cdr2': {
                    'contacts_file': str(cdr2_contacts_file),
                    'heatmaps_dir': str(cdr1_dir / 'cdr2_heatmaps'),
                },
                'cdr3': {
                    'contacts_file': str(cdr3_contacts_file),
                    'heatmaps_dir': str(cdr1_dir / 'cdr3_heatmaps'),
                },
                'non_cdr': {
                    'contacts_file': str(non_cdr_contacts_file),
                    'heatmaps_dir': str(cdr1_dir / 'non_cdr_heatmaps'),
                },
            }

            result = ContactHeatmapNode().execute(context)
            outputs = result.results['contact_heatmaps']

            self.assertTrue(Path(outputs['global']['peptide_tcr']['alpha']['image_file']).exists())
            self.assertIsNone(outputs['global']['peptide_tcr']['beta']['image_file'])
            self.assertTrue(Path(outputs['global']['hla_tcr']['beta']['image_file']).exists())
            self.assertIsNone(outputs['global']['hla_tcr']['alpha']['image_file'])
            self.assertTrue(Path(outputs['global']['groove_tcr']['beta']['image_file']).exists())
            self.assertIsNone(outputs['global']['groove_tcr']['alpha']['image_file'])
            self.assertTrue(Path(outputs['cdr1']['peptide_tcr']['alpha']['image_file']).exists())
            self.assertIsNone(outputs['cdr1']['hla_tcr']['alpha']['image_file'])
            self.assertIsNone(outputs['cdr1']['groove_tcr']['alpha']['image_file'])
            self.assertIsNone(outputs['cdr2']['peptide_tcr']['alpha']['image_file'])


if __name__ == '__main__':
    unittest.main()
