
db_testhead = ['SADI 0.02 C1 C2 C1 C3 C1 C4', 'SADI 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
               'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4 ',
               'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7 ',
               'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1 ',
               'SIMU O1 > F9', 'RIGU O1 > F9']
dbtest = {'ALO1': {'atoms': [['AL1', 'AL', '9.463', '-3.351', '3.397'],
                    ['O1', 'O', '8.422', '-2.079', '4.093'],
                    ['C1', 'C', '7.600', '-1.044', '4.188'],
                    ['C2', 'C', '6.402', '-1.006', '3.088'],
                    ['F1', 'F', '6.808', '-0.868', '1.814'],
                    ['F2', 'F', '5.537', '0.018', '3.271'],
                    ['F3', 'F', '5.634', '-2.113', '3.126'],
                    ['C3', 'C', '6.864', '-1.113', '5.645'],
                    ['F4', 'F', '7.745', '-1.050', '6.663'],
                    ['F5', 'F', '5.994', '-0.102', '5.866'],
                    ['F6', 'F', '6.133', '-2.224', '5.844'],
                    ['C4', 'C', '8.321', '0.415', '4.119'],
                    ['F7', 'F', '9.299', '0.599', '5.022'],
                    ['F8', 'F', '8.862', '0.674', '2.912'],
                    ['F9', 'F', '7.473', '1.448', '4.340'],
                    ['O2', 'O', '10.476', '-2.774', '2.046'],
                    ['C5', 'C', '11.284', '-2.995', '1.015'],
                    ['C6', 'C', '12.221', '-1.663', '0.860'],
                    ['F10', 'F', '13.116', '-1.734', '-0.152'],
                    ['F11', 'F', '11.498', '-0.555', '0.606'],
                    ['F12', 'F', '12.972', '-1.388', '1.941'],
                    ['C7', 'C', '10.554', '-3.196', '-0.430'],
                    ['F13', 'F', '9.842', '-4.339', '-0.502'],
                    ['F14', 'F', '11.424', '-3.271', '-1.464'],
                    ['F15', 'F', '9.723', '-2.199', '-0.779'],
                    ['C8', 'C', '12.317', '-4.246', '1.170'],
                    ['F16', 'F', '13.197', '-4.348', '0.147'],
                    ['F17', 'F', '13.085', '-4.139', '2.272'],
                    ['F18', 'F', '11.742', '-5.460', '1.217'],
                    ['O3', 'O', '8.442', '-4.660', '2.731'],
                    ['C9', 'C', '7.678', '-5.739', '2.656'],
                    ['C10', 'C', '6.896', '-5.722', '1.224'],
                    ['F19', 'F', '7.746', '-5.796', '0.181'],
                    ['F20', 'F', '6.140', '-4.628', '1.018'],
                    ['F21', 'F', '6.038', '-6.754', '1.061'],
                    ['C11', 'C', '8.439', '-7.177', '2.732'],
                    ['F22', 'F', '9.376', '-7.369', '1.786'],
                    ['F23', 'F', '9.047', '-7.386', '3.917'],
                    ['F24', 'F', '7.606', '-8.232', '2.581'],
                    ['C12', 'C', '6.524', '-5.782', '3.801'],
                    ['F25', 'F', '5.670', '-6.824', '3.672'],
                    ['F26', 'F', '5.736', '-4.689', '3.777'],
                    ['F27', 'F', '6.987', '-5.893', '5.059'],
                    ['O4', 'O', '10.460', '-3.990', '4.732'],
                    ['C14', 'C', '12.313', '-2.660', '5.662'],
                    ['C13', 'C', '11.237', '-3.873', '5.802'],
                    ['F28', 'F', '11.776', '-1.428', '5.645'],
                    ['F29', 'F', '13.057', '-2.764', '4.542'],
                    ['F30', 'F', '13.213', '-2.613', '6.672'],
                    ['C15', 'C', '10.516', '-3.668', '7.248'],
                    ['F31', 'F', '9.650', '-4.639', '7.590'],
                    ['F32', 'F', '9.843', '-2.503', '7.335'],
                    ['F33', 'F', '11.391', '-3.632', '8.281'],
                    ['C16', 'C', '12.124', '-5.241', '5.935'],
                    ['F34', 'F', '13.010', '-5.223', '6.957'],
                    ['F35', 'F', '11.358', '-6.327', '6.154'],
                    ['F36', 'F', '12.878', '-5.519', '4.856']],
          'comment': [['REM', 'Name:', 'AlOCCF334'],
                      ['REM', 'Produced', 'by', 'Grade', 'Web', 'Server',
                       'http://grade.globalphasing.org'],
                      ['REM', 'GEN:',
                       'Generated',
                       'by',
                       'GRADE',
                       '1.2.7',
                       '(February',
                       '19',
                       '2014)'],
                      ['REM', 'GEN:', 'from', 'mol2', 'file'],
                      ['REM',
                       'GEN:',
                       'using',
                       'quantum',
                       'mechanics',
                       'PM3'],
                      ['REM', 'grade-cif2shelx', 'output'],
                      ['REM',
                       'grade-cif2shelx',
                       'extracts',
                       'restraints',
                       'from',
                       'a',
                       'grade',
                       'CIF',
                       'file'],
                      ['REM', 'Version:', '0.0.5', '<Dec', '20', '2013>'],
                      ['REM', 'Total', 'charge', 'set', 'to', '-1']],
          'db': 'dsr-user-db',
          'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'],
          'head': [['DFIX', '1.785', '0.030', 'AL1', 'O1'],
                   ['DFIX', '1.784', '0.030', 'AL1', 'O2'],
                   ['DFIX', '1.788', '0.030', 'AL1', 'O3'],
                   ['DFIX', '1.785', '0.030', 'AL1', 'O4'],
                   ['DFIX', '1.325', '0.030', 'C1', 'O1'],
                   ['DFIX', '1.620', '0.030', 'C1', 'C2'],
                   ['DFIX', '1.628', '0.030', 'C1', 'C3'],
                   ['DFIX', '1.622', '0.030', 'C1', 'C4'],
                   ['DFIX', '1.344', '0.030', 'C2', 'F1'],
                   ['DFIX', '1.352', '0.030', 'C2', 'F2'],
                   ['DFIX', '1.347', '0.030', 'C2', 'F3'],
                   ['DFIX', '1.347', '0.030', 'C3', 'F4'],
                   ['DFIX', '1.351', '0.030', 'C3', 'F5'],
                   ['DFIX', '1.345', '0.030', 'C3', 'F6'],
                   ['DFIX', '1.344', '0.030', 'C4', 'F7'],
                   ['DFIX', '1.348', '0.030', 'C4', 'F8'],
                   ['DFIX', '1.352', '0.030', 'C4', 'F9'],
                   ['DFIX', '1.329', '0.030', 'C5', 'O2'],
                   ['DFIX', '1.630', '0.030', 'C5', 'C6'],
                   ['DFIX', '1.624', '0.030', 'C5', 'C7'],
                   ['DFIX', '1.622', '0.030', 'C5', 'C8'],
                   ['DFIX', '1.351', '0.030', 'C6', 'F10'],
                   ['DFIX', '1.347', '0.030', 'C6', 'F11'],
                   ['DFIX', '1.345', '0.030', 'C6', 'F12'],
                   ['DFIX', '1.348', '0.030', 'C7', 'F13'],
                   ['DFIX', '1.352', '0.030', 'C7', 'F14'],
                   ['DFIX', '1.344', '0.030', 'C7', 'F15'],
                   ['DFIX', '1.352', '0.030', 'C8', 'F16'],
                   ['DFIX', '1.347', '0.030', 'C8', 'F17'],
                   ['DFIX', '1.344', '0.030', 'C8', 'F18'],
                   ['DFIX', '1.324', '0.030', 'C9', 'O3'],
                   ['DFIX', '1.626', '0.030', 'C10', 'C9'],
                   ['DFIX', '1.622', '0.030', 'C11', 'C9'],
                   ['DFIX', '1.620', '0.030', 'C12', 'C9'],
                   ['DFIX', '1.347', '0.030', 'C10', 'F19'],
                   ['DFIX', '1.345', '0.030', 'C10', 'F20'],
                   ['DFIX', '1.351', '0.030', 'C10', 'F21'],
                   ['DFIX', '1.345', '0.030', 'C11', 'F22'],
                   ['DFIX', '1.348', '0.030', 'C11', 'F23'],
                   ['DFIX', '1.352', '0.030', 'C11', 'F24'],
                   ['DFIX', '1.352', '0.030', 'C12', 'F25'],
                   ['DFIX', '1.347', '0.030', 'C12', 'F26'],
                   ['DFIX', '1.345', '0.030', 'C12', 'F27'],
                   ['DFIX', '1.327', '0.030', 'C13', 'O4'],
                   ['DFIX', '1.622', '0.030', 'C13', 'C14'],
                   ['DFIX', '1.344', '0.030', 'C14', 'F28'],
                   ['DFIX', '1.348', '0.030', 'C14', 'F29'],
                   ['DFIX', '1.352', '0.030', 'C14', 'F30'],
                   ['DFIX', '1.622', '0.030', 'C13', 'C15'],
                   ['DFIX', '1.629', '0.030', 'C13', 'C16'],
                   ['DFIX', '1.345', '0.030', 'C15', 'F31'],
                   ['DFIX', '1.348', '0.030', 'C15', 'F32'],
                   ['DFIX', '1.352', '0.030', 'C15', 'F33'],
                   ['DFIX', '1.351', '0.030', 'C16', 'F34'],
                   ['DFIX', '1.347', '0.030', 'C16', 'F35'],
                   ['DFIX', '1.345', '0.030', 'C16', 'F36'],
                   ['DANG', '2.981', '0.062', 'O1', 'O2'],
                   ['DANG', '2.918', '0.064', 'O1', 'O3'],
                   ['DANG', '2.856', '0.066', 'O2', 'O3'],
                   ['DANG', '2.866', '0.065', 'O1', 'O4'],
                   ['DANG', '2.948', '0.063', 'O2', 'O4'],
                   ['DANG', '2.920', '0.064', 'O3', 'O4'],
                   ['DANG', '3.069', '0.044', 'AL1', 'C1'],
                   ['DANG', '2.502', '0.054', 'C2', 'O1'],
                   ['DANG', '2.408', '0.056', 'C3', 'O1'],
                   ['DANG', '2.577', '0.062', 'C2', 'C3'],
                   ['DANG', '2.501', '0.054', 'C4', 'O1'],
                   ['DANG', '2.577', '0.062', 'C2', 'C4'],
                   ['DANG', '2.578', '0.062', 'C3', 'C4'],
                   ['DANG', '2.506', '0.055', 'C1', 'F1'],
                   ['DANG', '2.478', '0.055', 'C1', 'F2'],
                   ['DANG', '2.127', '0.055', 'F1', 'F2'],
                   ['DANG', '2.472', '0.055', 'C1', 'F3'],
                   ['DANG', '2.159', '0.054', 'F1', 'F3'],
                   ['DANG', '2.138', '0.055', 'F2', 'F3'],
                   ['DANG', '2.474', '0.056', 'C1', 'F4'],
                   ['DANG', '2.492', '0.055', 'C1', 'F5'],
                   ['DANG', '2.145', '0.055', 'F4', 'F5'],
                   ['DANG', '2.505', '0.055', 'C1', 'F6'],
                   ['DANG', '2.158', '0.054', 'F4', 'F6'],
                   ['DANG', '2.129', '0.055', 'F5', 'F6'],
                   ['DANG', '2.503', '0.055', 'C1', 'F7'],
                   ['DANG', '2.480', '0.055', 'C1', 'F8'],
                   ['DANG', '2.160', '0.054', 'F7', 'F8'],
                   ['DANG', '2.482', '0.055', 'C1', 'F9'],
                   ['DANG', '2.126', '0.055', 'F7', 'F9'],
                   ['DANG', '2.138', '0.055', 'F8', 'F9'],
                   ['DANG', '3.019', '0.046', 'AL1', 'C5'],
                   ['DANG', '2.391', '0.057', 'C6', 'O2'],
                   ['DANG', '2.519', '0.054', 'C7', 'O2'],
                   ['DANG', '2.582', '0.062', 'C6', 'C7'],
                   ['DANG', '2.519', '0.054', 'C8', 'O2'],
                   ['DANG', '2.575', '0.062', 'C6', 'C8'],
                   ['DANG', '2.573', '0.062', 'C7', 'C8'],
                   ['DANG', '2.496', '0.055', 'C5', 'F10'],
                   ['DANG', '2.477', '0.056', 'C5', 'F11'],
                   ['DANG', '2.140', '0.055', 'F10', 'F11'],
                   ['DANG', '2.505', '0.055', 'C5', 'F12'],
                   ['DANG', '2.127', '0.055', 'F10', 'F12'],
                   ['DANG', '2.158', '0.054', 'F11', 'F12'],
                   ['DANG', '2.482', '0.055', 'C5', 'F13'],
                   ['DANG', '2.479', '0.055', 'C5', 'F14'],
                   ['DANG', '2.138', '0.055', 'F13', 'F14'],
                   ['DANG', '2.505', '0.055', 'C5', 'F15'],
                   ['DANG', '2.164', '0.054', 'F13', 'F15'],
                   ['DANG', '2.126', '0.055', 'F14', 'F15'],
                   ['DANG', '2.479', '0.055', 'C5', 'F16'],
                   ['DANG', '2.471', '0.055', 'C5', 'F17'],
                   ['DANG', '2.138', '0.055', 'F16', 'F17'],
                   ['DANG', '2.511', '0.054', 'C5', 'F18'],
                   ['DANG', '2.123', '0.055', 'F16', 'F18'],
                   ['DANG', '2.162', '0.054', 'F17', 'F18'],
                   ['DANG', '3.071', '0.044', 'AL1', 'C9'],
                   ['DANG', '2.411', '0.056', 'C10', 'O3'],
                   ['DANG', '2.521', '0.054', 'C11', 'O3'],
                   ['DANG', '2.579', '0.062', 'C10', 'C11'],
                   ['DANG', '2.472', '0.055', 'C12', 'O3'],
                   ['DANG', '2.582', '0.062', 'C10', 'C12'],
                   ['DANG', '2.574', '0.062', 'C11', 'C12'],
                   ['DANG', '2.471', '0.056', 'C9', 'F19'],
                   ['DANG', '2.503', '0.055', 'C9', 'F20'],
                   ['DANG', '2.156', '0.054', 'F19', 'F20'],
                   ['DANG', '2.490', '0.055', 'C9', 'F21'],
                   ['DANG', '2.148', '0.054', 'F19', 'F21'],
                   ['DANG', '2.129', '0.055', 'F20', 'F21'],
                   ['DANG', '2.507', '0.055', 'C9', 'F22'],
                   ['DANG', '2.480', '0.055', 'C9', 'F23'],
                   ['DANG', '2.159', '0.054', 'F22', 'F23'],
                   ['DANG', '2.479', '0.055', 'C9', 'F24'],
                   ['DANG', '2.125', '0.055', 'F22', 'F24'],
                   ['DANG', '2.139', '0.055', 'F23', 'F24'],
                   ['DANG', '2.483', '0.055', 'C9', 'F25'],
                   ['DANG', '2.470', '0.055', 'C9', 'F26'],
                   ['DANG', '2.138', '0.055', 'F25', 'F26'],
                   ['DANG', '2.502', '0.055', 'C9', 'F27'],
                   ['DANG', '2.130', '0.055', 'F25', 'F27'],
                   ['DANG', '2.161', '0.054', 'F26', 'F27'],
                   ['DANG', '3.034', '0.045', 'AL1', 'C13'],
                   ['DANG', '2.506', '0.055', 'C13', 'F28'],
                   ['DANG', '2.471', '0.055', 'C13', 'F29'],
                   ['DANG', '2.157', '0.054', 'F28', 'F29'],
                   ['DANG', '2.484', '0.055', 'C13', 'F30'],
                   ['DANG', '2.129', '0.055', 'F28', 'F30'],
                   ['DANG', '2.141', '0.055', 'F29', 'F30'],
                   ['DANG', '2.469', '0.055', 'C14', 'O4'],
                   ['DANG', '2.541', '0.053', 'C15', 'O4'],
                   ['DANG', '2.579', '0.062', 'C14', 'C15'],
                   ['DANG', '2.410', '0.056', 'C16', 'O4'],
                   ['DANG', '2.576', '0.062', 'C14', 'C16'],
                   ['DANG', '2.577', '0.062', 'C15', 'C16'],
                   ['DANG', '2.507', '0.055', 'C13', 'F31'],
                   ['DANG', '2.480', '0.055', 'C13', 'F32'],
                   ['DANG', '2.163', '0.054', 'F31', 'F32'],
                   ['DANG', '2.477', '0.055', 'C13', 'F33'],
                   ['DANG', '2.125', '0.055', 'F31', 'F33'],
                   ['DANG', '2.138', '0.055', 'F32', 'F33'],
                   ['DANG', '2.493', '0.055', 'C13', 'F34'],
                   ['DANG', '2.477', '0.056', 'C13', 'F35'],
                   ['DANG', '2.143', '0.055', 'F34', 'F35'],
                   ['DANG', '2.506', '0.055', 'C13', 'F36'],
                   ['DANG', '2.127', '0.055', 'F34', 'F36'],
                   ['DANG', '2.158', '0.054', 'F35', 'F36'],
                   ['SIMU', 'AL1', '>', 'F36'],
                   ['RIGU', 'AL1', '>', 'F36']],
          'line': None,
          'name': 'ALO1',
          'resi': 'ALO1'}}
