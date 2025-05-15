#----------------change file path-------------------#
import sys
import os

# Go up two levels and then down into 'target_module'
src_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..'))

test_dir = src_dir + '/tests'

# Add that directory to sys.path
sys.path.insert(0, src_dir)
sys.path.insert(0, test_dir)

# check directory is right
print(src_dir)
#-------------------------------------------------#

from datetime import datetime
import pandas as pd
import metalsilicate as ms

pd.set_option('display.max_rows', None, 'display.max_columns', None)

myfile = ms.ExcelFile(test_dir + '/testIn_data_three-samples.xlsx')

#Settings for Saved File Name
todays_date = datetime.today().strftime('%Y-%m-%d')
filename = test_dir + '/TEST_out_' + '.xlsx'

#Process data
processed = myfile.process(pressure="P (GPa)", temperature="T_C", filename=filename)

#-----------------DO THE TEST---------------------#
import unittest
import pathlib

#####
def print_msg_box(msg, indent=1, width=None, title=None):
    """Print message-box with optional title."""
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print("\n")
    print(box)
#####

TEST_FILE = pathlib.Path(__file__).parent.joinpath("TEST_out_.xlsx")
TRUTH_FILE = pathlib.Path(__file__).parent.joinpath("truthOut_three-samples.xlsx")

class TestOutputValuesBatch(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg):
        """
        Creates a new type of unittest to assert that pd DataFrames
        are equal, inheriting from pandas testing routine.
        If not, prints the diff between them.
        """
        try:
            pd._testing.assert_frame_equal(a, b, check_dtype=False)
        except AssertionError as e:
            print_msg_box("DataFrames differ! Showing difference:")

            # Show where they differ
            try:
                comparison = self.test_df.compare(
                self.truth_df,
                keep_shape=False, # does not display colums that are equal
                keep_equal=False,
                result_names=("Output", "Truth"))
                print(comparison)
            except Exception as compare_error:
                print("Error comparing dataframes:", compare_error)

            raise self.failureException(f"{msg}\n{str(e)}") from e
    
    def setUp(self):
        # Add assertDataframeEqual to unittest test cases
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

        # import test file to pandas df
        self.test_df = pd.read_excel(TEST_FILE, sheet_name="Thermo_Params")

        # import truth file to compare test to as pandas df
        self.truth_df = pd.read_excel(TRUTH_FILE)
    
    def test_ImportExcel(self):
        print_msg_box("TestOutputValuesBatch")
        self.assertEqual(self.test_df, self.truth_df, 
                         'DataFrames are different')

#-------------------------------------------------#

if __name__ == '__main__':
    unittest.main()
