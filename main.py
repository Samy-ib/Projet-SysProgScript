import gui
import sys
from PyQt5.QtWidgets import QFileDialog
from fonc import *
# output_result = ""

# input_list = []

# def print_output(liste):
#     #FOnction qui servira a afficher les rsulatts



class fenetre(gui.Ui_Form):
    
    def setupUi(self, Form):
        super().setupUi(Form)  

        self.output_result = ""
        self.input_list = ""
        self.input_seq = ""

        self.button_import.clicked.connect(self.importer)
        self.button_generate.clicked.connect(self.generate)
        self.button_valider.clicked.connect(self.valider_in)
        self.comboBox.currentIndexChanged.connect(self.combo_lineEdit)
        self.in_textEdit.textChanged.connect(self.textChange)


        self.button_validite.clicked.connect(self.check_validite)
        self.button_frequence.clicked.connect(self.frequences)
        self.button_adnArn.clicked.connect(self.toArn)
        self.button_arnProteine.clicked.connect(self.toProteine)
        self.button_compInver.clicked.connect(self.compInv)
        self.button_taux_gc.clicked.connect(self.tauxGc)
        self.button_freqCodon.clicked.connect(self.freqCodon)
        self.button_masseProteique.clicked.connect(self.masseProteique)
        self.button_epissage.clicked.connect(self.epissage)
        self.button_assemblage.clicked.connect(self.assemblage)
        self.button_sauvegarder.clicked.connect(self.save)

        ############# INPUT HANDEL BUTTONS #############
    
    def combo_lineEdit(self):
        self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#f7b12c;">En attente de validation</span></p></body></html>')
        if len(self.comboBox.currentText())>0 :
            self.in_textEdit.setText(str(self.input_list[int(self.comboBox.currentText()[0])-1][1]))

    def importer(self):
        filename =QFileDialog.getOpenFileName(MainWindow,"QFileDialog.getOpenFileName()", "","Fichier Fasta (*.fas)")
        if(filename[0]):
            self.in_textEdit.setText(str(self.parse_fasta(filename[0])))
            self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#f7b12c;;">En attente de validation</span></p></body></html>')
            self.comboBox.clear()
            for i in range(len(self.input_list)):
                val = str(i+1) + "-" + str(self.input_list[i][0]) 
                self.comboBox.addItem(val)

    def parse_fasta(self, filename):
        self.input_list = [[rec.id, rec.seq] for rec in SeqIO.parse(filename, "fasta")]
        return self.input_list[0][1]

    def generate(self):
        self.in_textEdit.setText(generate())
        self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#f7b12c;;">En attente de validation</span></p></body></html>')
    
    def valider_in(self):
        self.input_seq = (self.in_textEdit.toPlainText()).upper().replace("\n","").replace(" ","")

        if(not self.in_textEdit.toPlainText()):
            self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#AF3D35;;">Aucune sequence !</span></p></body></html>')
        elif not valide(self.input_seq):
            self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#AF3D35;;">Sequence non valide !</span></p></body></html>')
        else:
            self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#f7b12c;">Sequence validee !</span></p></body></html>')
            self.button_validite.setEnabled(True)
            self.button_frequence.setEnabled(True)
            self.button_adnArn.setEnabled(True)
            self.button_arnProteine.setEnabled(True)
            self.button_compInver.setEnabled(True)
            self.button_taux_gc.setEnabled(True)
            self.button_freqCodon.setEnabled(True)
            self.button_masseProteique.setEnabled(True)
            self.button_epissage.setEnabled(True)
            self.button_assemblage.setEnabled(True)

        self.out_textEdit.clear()
        

        self.output_result = ""

    def textChange(self):
        self.label_valider.setText('<html><head/><body><p align="center"><span style=" font-size:18pt; color:#f7b12c;">En attente de validation</span></p></body></html>')

    ############# OUTPUT HANDEL BUTTONS #############

    def check_validite(self):
        if valide(self.input_seq):
            out="La sequence est une sequence valide !"
        else:
            out="La sequence n'est pas valide !"  
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def frequences(self):
        out = calcul_freq(self.input_seq)
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def toArn(self):
        arn = adn_to_arn(self.input_seq)
        out = "ARN : " + arn
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)
    
    def toProteine(self):
        prot = adn_to_proteine(self.input_seq)
        out = "Proteine : " + prot
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def compInv(self):
        ci = comp_inv(list(self.input_seq))
        out = "Complement inverse : " + ci
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def tauxGc(self):
        tgc = taux_gc(self.input_seq)
        out = "Taux de GC : " + str("%.2f" % tgc) + "%"
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)
    
    def freqCodon(self):
        out = "Frequence des codons :"+ "\n" + freq_codon(self.input_seq)
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def masseProteique(self):
        masse = masse_proteique(self.input_seq)
        out = "Masse proteique : " + str(masse) + "Da"
        if out in self.output_result:
            self.output_result = self.output_result.replace(out + "\n\n","")
        else:
            self.output_result += out + "\n\n"
        self.out_textEdit.setText(self.output_result)

    def epissage(self):
        if (self.lineEdit.text().upper() in adn_to_arn(self.input_seq)) and (self.lineEdit_2.text().upper() in adn_to_arn(self.input_seq)) and self.lineEdit.text() and self.lineEdit_2.text() :
            epi = epissage_adn(adn_to_arn(self.input_seq), self.lineEdit.text().upper(), self.lineEdit_2.text().upper())
            self.label_introns.setText('<html><head/><body><p>Introns: <span style=" color:#f7b12c;">Correct !</span></p></body></html>')
            out = "Epissage ADN : " + epi
            if out in self.output_result:
                self.output_result = self.output_result.replace(out + "\n\n","")
            else:
                self.output_result += out + "\n\n"
            self.out_textEdit.setText(self.output_result)
        else:
            self.label_introns.setText('<html><head/><body><p>Introns: <span style=" color:#e82b1f;">Erreur ! </span></p></body></html>')

    def assemblage(self):
        try:
            taille = int(self.lineEdit_4.text())
            self.label_4.setText('<html><head/><body><p>Taille:</p></body></html>')
            seq = assem(self.input_seq, taille)
            out = "Assemblage : " + seq
            if out in self.output_result:
                self.output_result = self.output_result.replace(out + "\n\n","")
            else:
                self.output_result += out + "\n\n"
            self.out_textEdit.setText(self.output_result)
        except ValueError:
            self.label_4.setText('<html><head/><body><p>Taille: <span style=" color:#e8291f;">Erreur</span></p></body></html>')
        
    def save(self):
        if(self.lineEdit_3.text()):
            fichier = open("./saves/"+self.lineEdit_3.text()+".txt", "w")
            fichier.write(self.output_result)
            fichier.close()
if __name__ == '__main__':
    app = gui.QtWidgets.QApplication(sys.argv)
    MainWindow = gui.QtWidgets.QMainWindow()
    MainWindow.setWindowFlags(gui.QtCore.Qt.FramelessWindowHint)
    ui = fenetre()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
