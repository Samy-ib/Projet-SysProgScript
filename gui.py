# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'untitled.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1024, 600)
        Form.setMinimumSize(QtCore.QSize(1024, 600))
        Form.setMaximumSize(QtCore.QSize(1024, 600))
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        font.setKerning(False)
        Form.setFont(font)
        Form.setStyleSheet("border-radius: 4px;\n"
"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0 rgba(27, 55, 68, 255), stop:1 rgba(0, 171, 164, 255));\n"
"")
        self.button_generate = QtWidgets.QPushButton(Form)
        self.button_generate.setGeometry(QtCore.QRect(20, 330, 141, 31))
        self.button_generate.setStyleSheet("background-color: rgba(7, 163, 148,0.5);\n"
"color: rgba(0, 0, 0,0.8)")
        self.button_generate.setObjectName("button_generate")
        self.button_import = QtWidgets.QPushButton(Form)
        self.button_import.setGeometry(QtCore.QRect(200, 330, 141, 31))
        self.button_import.setStyleSheet("background-color: rgba(7, 163, 148,0.5);color: rgba(0, 0, 0,0.8)")
        self.button_import.setObjectName("button_import")
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(10, 30, 121, 20))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.label.setFont(font)
        self.label.setStyleSheet("background: rgba(0,0,0,0.01);\n"
"color: rgba(166, 212, 194,0.5);")
        self.label.setObjectName("label")
        self.button_validite = QtWidgets.QPushButton(Form)
        self.button_validite.setEnabled(False)
        self.button_validite.setGeometry(QtCore.QRect(850, 60, 141, 31))
        self.button_validite.setStyleSheet("background-color: rgba(7, 163, 148,0.5);color: rgba(0, 0, 0,0.8);")
        self.button_validite.setObjectName("button_validite")
        self.button_frequence = QtWidgets.QPushButton(Form)
        self.button_frequence.setEnabled(False)
        self.button_frequence.setGeometry(QtCore.QRect(850, 100, 141, 31))
        self.button_frequence.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_frequence.setObjectName("button_frequence")
        self.button_adnArn = QtWidgets.QPushButton(Form)
        self.button_adnArn.setEnabled(False)
        self.button_adnArn.setGeometry(QtCore.QRect(850, 140, 141, 31))
        self.button_adnArn.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_adnArn.setObjectName("button_adnArn")
        self.button_arnProteine = QtWidgets.QPushButton(Form)
        self.button_arnProteine.setEnabled(False)
        self.button_arnProteine.setGeometry(QtCore.QRect(850, 180, 141, 31))
        self.button_arnProteine.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_arnProteine.setObjectName("button_arnProteine")
        self.button_compInver = QtWidgets.QPushButton(Form)
        self.button_compInver.setEnabled(False)
        self.button_compInver.setGeometry(QtCore.QRect(850, 220, 141, 31))
        self.button_compInver.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_compInver.setObjectName("button_compInver")
        self.button_taux_gc = QtWidgets.QPushButton(Form)
        self.button_taux_gc.setEnabled(False)
        self.button_taux_gc.setGeometry(QtCore.QRect(850, 260, 141, 31))
        self.button_taux_gc.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_taux_gc.setObjectName("button_taux_gc")
        self.button_freqCodon = QtWidgets.QPushButton(Form)
        self.button_freqCodon.setEnabled(False)
        self.button_freqCodon.setGeometry(QtCore.QRect(850, 300, 141, 31))
        self.button_freqCodon.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_freqCodon.setObjectName("button_freqCodon")
        self.button_masseProteique = QtWidgets.QPushButton(Form)
        self.button_masseProteique.setEnabled(False)
        self.button_masseProteique.setGeometry(QtCore.QRect(850, 340, 141, 31))
        self.button_masseProteique.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_masseProteique.setObjectName("button_masseProteique")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(460, 30, 121, 20))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.label_2.setFont(font)
        self.label_2.setStyleSheet("background: rgba(0,0,0,0);\n"
"color: rgba(166, 212, 194,0.5);")
        self.label_2.setObjectName("label_2")
        self.button_epissage = QtWidgets.QPushButton(Form)
        self.button_epissage.setEnabled(False)
        self.button_epissage.setGeometry(QtCore.QRect(850, 460, 141, 31))
        self.button_epissage.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_epissage.setObjectName("button_epissage")
        self.button_assemblage = QtWidgets.QPushButton(Form)
        self.button_assemblage.setEnabled(False)
        self.button_assemblage.setGeometry(QtCore.QRect(850, 510, 141, 31))
        self.button_assemblage.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_assemblage.setObjectName("button_assemblage")
        self.button_sauvegarder = QtWidgets.QPushButton(Form)
        self.button_sauvegarder.setGeometry(QtCore.QRect(690, 560, 141, 31))
        self.button_sauvegarder.setStyleSheet("color: rgba(0, 0, 0,0.8);background-color: rgba(7, 163, 148,0.5);")
        self.button_sauvegarder.setObjectName("button_sauvegarder")
        self.validite_button_2 = QtWidgets.QPushButton(Form)
        self.validite_button_2.setGeometry(QtCore.QRect(980, 10, 31, 31))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        self.validite_button_2.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Noto Mono")
        font.setPointSize(17)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.validite_button_2.setFont(font)
        self.validite_button_2.setStyleSheet("background-color: rgba(232, 43, 31,0);\n"
"color: rgb(232, 43, 31);")
        self.validite_button_2.setObjectName("validite_button_2")
        self.validite_button_3 = QtWidgets.QPushButton(Form)
        self.validite_button_3.setGeometry(QtCore.QRect(940, 10, 31, 31))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(232, 43, 31, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        self.validite_button_3.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Noto Mono")
        font.setPointSize(24)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.validite_button_3.setFont(font)
        self.validite_button_3.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.validite_button_3.setAutoFillBackground(False)
        self.validite_button_3.setStyleSheet("background-color: rgba(232, 43, 31,0);\n"
"color: rgb(232, 43, 31);")
        self.validite_button_3.setObjectName("validite_button_3")
        self.in_textEdit = QtWidgets.QTextEdit(Form)
        self.in_textEdit.setGeometry(QtCore.QRect(10, 60, 341, 261))
        self.in_textEdit.setStyleSheet("border-radius:4px;\n"
"background-color: rgba(0, 151, 167,0.5);\n"
"color: rgb(166, 212, 194);")
        self.in_textEdit.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.in_textEdit.setReadOnly(False)
        self.in_textEdit.setCursorWidth(0)
        self.in_textEdit.setObjectName("in_textEdit")
        self.out_textEdit = QtWidgets.QTextEdit(Form)
        self.out_textEdit.setGeometry(QtCore.QRect(460, 60, 371, 481))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.NoBrush)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.NoBrush)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(166, 212, 194))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.NoBrush)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 151, 167, 127))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        self.out_textEdit.setPalette(palette)
        self.out_textEdit.setStyleSheet("border-radius:4px; background-color: rgba(0, 151, 167,0.5);\n"
"color: rgb(166, 212, 194);")
        self.out_textEdit.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.out_textEdit.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.out_textEdit.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.out_textEdit.setReadOnly(True)
        self.out_textEdit.setObjectName("out_textEdit")
        self.button_valider = QtWidgets.QPushButton(Form)
        self.button_valider.setGeometry(QtCore.QRect(110, 380, 141, 31))
        self.button_valider.setStyleSheet("background-color: rgba(7, 163, 148,0.5);color: rgba(0, 0, 0,0.8)")
        self.button_valider.setCheckable(False)
        self.button_valider.setObjectName("button_valider")
        self.label_3 = QtWidgets.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(10, 490, 361, 81))
        self.label_3.setStyleSheet("background: rgba(0,0,0,0.01);\n"
"color: rgba(166, 212, 194,1);")
        self.label_3.setObjectName("label_3")
        self.label_valider = QtWidgets.QLabel(Form)
        self.label_valider.setGeometry(QtCore.QRect(30, 420, 311, 51))
        self.label_valider.setStyleSheet("background: rgba(0,0,0,0.01);")
        self.label_valider.setObjectName("label_valider")
        self.comboBox = QtWidgets.QComboBox(Form)
        self.comboBox.setGeometry(QtCore.QRect(360, 90, 83, 29))
        self.comboBox.setStyleSheet("background-color: rgba(0, 151, 167,0.5);")
        self.comboBox.setObjectName("comboBox")
        self.label_introns = QtWidgets.QLabel(Form)
        self.label_introns.setGeometry(QtCore.QRect(850, 380, 171, 20))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.label_introns.setFont(font)
        self.label_introns.setStyleSheet("background: rgba(0,0,0,0);\n"
"color: rgba(166, 212, 194,0.5);")
        self.label_introns.setObjectName("label_introns")
        self.lineEdit = QtWidgets.QLineEdit(Form)
        self.lineEdit.setGeometry(QtCore.QRect(850, 400, 131, 21))
        self.lineEdit.setStyleSheet("border-radius:4px;\n"
"background-color: rgba(0, 151, 167,0.9);\n"
"color: rgb(166, 212, 194);")
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit_2 = QtWidgets.QLineEdit(Form)
        self.lineEdit_2.setGeometry(QtCore.QRect(850, 430, 131, 21))
        self.lineEdit_2.setStyleSheet("border-radius:4px;\n"
"background-color: rgba(0, 151, 167,0.9);\n"
"color: rgb(166, 212, 194);")
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.lineEdit_3 = QtWidgets.QLineEdit(Form)
        self.lineEdit_3.setGeometry(QtCore.QRect(460, 560, 221, 21))
        self.lineEdit_3.setStyleSheet("border-radius:4px;\n"
"background-color: rgba(0, 151, 167,0.9);\n"
"color: rgb(166, 212, 194);")
        self.lineEdit_3.setText("")
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.label_introns_2 = QtWidgets.QLabel(Form)
        self.label_introns_2.setGeometry(QtCore.QRect(460, 540, 271, 20))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.label_introns_2.setFont(font)
        self.label_introns_2.setStyleSheet("background: rgba(0,0,0,0);\n"
"color: rgba(166, 212, 194,0.5);")
        self.label_introns_2.setObjectName("label_introns_2")
        self.label_4 = QtWidgets.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(850, 490, 171, 20))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.label_4.setFont(font)
        self.label_4.setStyleSheet("background: rgba(0,0,0,0);\n"
"color: rgba(166, 212, 194,0.5);")
        self.label_4.setObjectName("label_4")
        self.lineEdit_4 = QtWidgets.QLineEdit(Form)
        self.lineEdit_4.setGeometry(QtCore.QRect(960, 490, 31, 21))
        self.lineEdit_4.setStyleSheet("border-radius:4px;\n"
"background-color: rgba(0, 151, 167,0.9);\n"
"color: rgb(166, 212, 194);")
        self.lineEdit_4.setMaxLength(3)
        self.lineEdit_4.setObjectName("lineEdit_4")

        self.retranslateUi(Form)
        self.validite_button_2.clicked.connect(Form.close)
        self.validite_button_3.clicked.connect(Form.lower)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.button_generate.setText(_translate("Form", "Genere aleatoirement"))
        self.button_import.setText(_translate("Form", "Importer un fichier"))
        self.label.setText(_translate("Form", "Votre sequence :"))
        self.button_validite.setText(_translate("Form", "Validite"))
        self.button_frequence.setText(_translate("Form", "Frequences"))
        self.button_adnArn.setText(_translate("Form", "ADN -> ARN"))
        self.button_arnProteine.setText(_translate("Form", "ARN -> Proteines"))
        self.button_compInver.setText(_translate("Form", "Complement inverse"))
        self.button_taux_gc.setText(_translate("Form", "Taux de GC"))
        self.button_freqCodon.setText(_translate("Form", "Frequences de codons"))
        self.button_masseProteique.setText(_translate("Form", "Masse proteique"))
        self.label_2.setText(_translate("Form", "Votre resultat :"))
        self.button_epissage.setText(_translate("Form", "Epissage ARN"))
        self.button_assemblage.setText(_translate("Form", "Assemblage"))
        self.button_sauvegarder.setText(_translate("Form", "Sauvegarder"))
        self.validite_button_2.setText(_translate("Form", "X"))
        self.validite_button_3.setText(_translate("Form", "-"))
        self.in_textEdit.setHtml(_translate("Form", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Inconsolata\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:456;\">ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG</span></p></body></html>"))
        self.out_textEdit.setHtml(_translate("Form", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Inconsolata\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:456;\"><br /></p></body></html>"))
        self.button_valider.setText(_translate("Form", "Valider"))
        self.label_3.setText(_translate("Form", "<html><head/><body><p>Veuillez soit generer une sequence ADN aleatoire </p><p>ou bien importer un fichier fasta.</p><p>Une fois ceci fait, cliquez sur les bouton &quot;valider&quot;.</p></body></html>"))
        self.label_valider.setText(_translate("Form", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; color:#f7b12c;;\">En attente de validation</span></p></body></html>"))
        self.label_introns.setText(_translate("Form", "<html><head/><body><p>Introns:</p></body></html>"))
        self.lineEdit.setText(_translate("Form", "AUCGGUCGAA"))
        self.lineEdit_2.setText(_translate("Form", "AUCGGUCGAGCGUGU"))
        self.label_introns_2.setText(_translate("Form", "<html><head/><body><p>Nom de votre fichier :</p></body></html>"))
        self.label_4.setText(_translate("Form", "<html><head/><body><p>Taille:</p></body></html>"))

