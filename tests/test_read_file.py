"""
TESTNI MODUL

Ucitavanje spremljene topologije iz datoteke i provjera

    1. Dio: Provjera opcih metoda definiranih u modulu grillage_model
    2. Dio: Provjera karakteristika ucitane topologije

Za provjeru određene metode ili karakteristike, u listi provjera ukloniti oznaku komentara # sa Ctrl + /
"""

from grillage.grillage_model import *
from timeit import default_timer as timer
np.set_printoptions(suppress=True)          # Suppress scientific output notation

start = timer()

hc_var = 1
filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
hc_variant = GrillageModelData(filename).read_file()        # Učitavanje topologije iz datoteke

matST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
FB_beam = FBBeamProperty(6, 1089, 10, matST24)

hc_variant.plating()[6].set_intercostal_stiffeners(4, FB_beam)    # Dodavanje interkostalnih ukrepa na sva neukrepljena polja zone 6
# hc_variant.plating()[7].set_intercostal_stiffeners(4, FB_beam)
# hc_variant.plating()[10].set_intercostal_stiffeners(4, FB_beam)
# hc_variant.plating()[11].set_intercostal_stiffeners(4, FB_beam)

# Pojedinačno dodavanje interkostalnih ukrepa - na treći elementarni panel druge zone oplate
# hc_variant.plating()[2].elementary_plate_panels[3].intercostal_stiffener_num = 1      # Jedna interkostalna ukrepa
# hc_variant.plating()[2].elementary_plate_panels[3].beam_prop = FB_beam                # beam property interkostala

# Ponovno spremanje ucitane topologije pod novim imenom "read_test.txt"
# GrillageModelData("read_test.txt").write_file(hc_variant)


def KoordinateCvorova():
    # Koordinate sjecista proizvoljno odabranih tocaka
    # Eksplicitno zadani cvorovi krajeva dva okomita nosaca (id, x, y, z) za provjeru:

    # Cvorovi koji definiraju prvi nosac:
    node1 = ModelNode(1, 0, 10, 0)
    node2 = ModelNode(2, 18, 10, 0)

    # Cvorovi koji definiraju drugi nosac:
    node3 = ModelNode(3, 12, 0, 0)
    node4 = ModelNode(4, 12, 20, 0)

    # Vracanje koordinata kroz metodu coords
    n1 = node1.coords
    n2 = node2.coords
    n3 = node3.coords
    n4 = node4.coords

    print("Uneseni cvorovi prvog nosaca:", n1, ",", n2)
    print("Uneseni cvorovi drugog nosaca:", n3, ",", n4)
    print("Koordinate sjecista nosaca: ", Grillage.get_intersection(n1, n2, n3, n4))


def SvojstvaTprofila(hw, tw, bf, tf, bp, tp):
    #   Provjera karakteristika proizvoljno unesenog T profila

    # Eksplicitno zadavanje karakteristika T profila:
    #    (id, hw, tw, bf, tf, mat)

    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    beam = TBeamProperty(1, hw, tw, bf, tf, ST24)
    tc1 = CorrosionAddition(1, 2)

    # beam1.hw = 69         # Izmjena visine struka
    print(" Visina struka (gross),                          hw = ", beam.hw, "mm")
    print(" Debljina struka (gross),                        tw = ", beam.tw, "mm")
    print(" Sirina prirubnice (gross),                      bf = ", beam.bf, "mm")
    print(" Debljina prirubnice (gross),                    tf = ", beam.tf, "mm")
    print(" Sirina sunosivog oplocenja,                     bp = ", bp, "mm")
    print(" Debljina sunosivog oplocenja (gross),           tp = ", tp, "mm")
    print(" Korozijski dodatak,                             tc = ", tc1.tc, "mm")
    print(" Visina struka (net),                        hw_net = ", beam.hw_net(tc1, tp))
    print(" Debljina struka (net),                      tw_net = ", beam.tw_net(tc1), "mm")
    print(" Sirina prirubnice (net),                    bf_net = ", beam.bf_net(tc1), "mm")
    print(" Debljina prirubnice (net),                  tf_net = ", beam.tf_net(tc1), "mm")
    print(" Debljina sunosivog oplocenja (net),         tp_net = ", PlateProperty.tp_net(tc1, tp), "mm")
    print(" Granica razvlacenja,                           Reh = ", beam.mat.Reh, "N/mm2")
    print(" Smicna povrsina (net),                        A_sh = ", beam.getShArea_T(tp, tc1), "cm2")
    print(" Povrsina samog T profila (gross),              A_T = ", "{:.2f}".format(beam.getArea), "cm2")
    print(" Povrsina nosaca sa sunosivom sirinom (net),      A = ", "{:.2f}".format(beam.getArea_I(bp, tp, tc1)), "cm2")
    print(" Moment inercije (net),                          Iy = ", "{:.2f}".format(beam.get_Iy_I(bp, tp, tc1)), "cm4")
    print(" Moment otpora (net),                          Wmin = ", "{:.2f}".format(beam.getWmin(bp, tp, tc1)), "cm3")
    print(" Sektorski moment inercije - warping (net),      Iw = ", "{:.2f}".format(beam.get_Iw(tc1)), "cm6")
    print(" Polarni moment inercije (net),                  Ip = ", "{:.2f}".format(beam.get_Ip(tc1)), "cm4")
    print(" St Venant-ov moment inercije (net),             It = ", "{:.2f}".format(beam.get_It(tc1)), "cm4")


def SvojstvaHPprofila(h_HP, t_HP, bp, tp):
    #   Provjera karakteristika proizvoljno unesenog HP profila

    # Eksplicitno zadavanje karakteristika HP profila:
    AH36 = MaterialProperty(1, 210000, 0.3, 7850, 355, "AH36")
    hp1 = BulbBeamProperty(1, h_HP, t_HP, AH36)
    tc1 = CorrosionAddition(1, 2)

    print(" Visina HP profila,                               hw_HP = ", hp1.hw_HP, "mm")
    print(" Debljina struka HP profila,                      tw_HP = ", hp1.tw_HP, "mm")
    print(" Visina struka ekvivalentnog L (gross),              hw = ", "{:.2f}".format(hp1.hw_ekv), "mm")
    print(" Debljina struka ekvivalentnog L (gross),            tw = ", hp1.tw_ekv, "mm")
    print(" Sirina prirubnice ekvivalentnog L (gross),          bf = ", "{:.2f}".format(hp1.bf_ekv), "mm")
    print(" Debljina prirubnice ekvivalentnog L (gross),        tf = ", "{:.2f}".format(hp1.tf_ekv), "mm")
    print(" Sirina sunosivog oplocenja,                         bp = ", bp, "mm")
    print(" Debljina sunosivog oplocenja  (gross),              tp = ", tp, "mm")
    print(" Korozijski dodatak,                                 tc = ", tc1.tc, "mm")
    print(" Granica razvlacenja,                               Reh = ", hp1.mat.Reh, "N/mm2")
    print(" Gustoca materijala,                                 ro = ", hp1.mat.ro, "kg/m3")
    print(" Smicna povrsina (net),                            A_sh = ", "{:.2f}".format(hp1.getShArea_HP(tc1)), "cm2")
    print(" Povrsina samog HP profila (gross),                A_HP = ", "{:.2f}".format(hp1.getArea), "cm2")
    print(" Povrsina HP profila sa sunosivom sirinom (net),      A = ", "{:.2f}".format(hp1.getArea_I(bp, tp, tc1)), "cm2")
    print(" Težište HP profila od oplate (net),               z_NA = ", "{:.2f}".format(hp1.get_z_na_I(bp, tp, tc1)), "mm")
    print(" Moment inercije (net),                              Iy = ", "{:.2f}".format(hp1.get_Iy_I(bp, tp, tc1)), "cm4")
    print(" Moment otpora (net),                              Wmin = ", "{:.2f}".format(hp1.getWmin(bp, tp, tc1)), "cm3")
    print(" Sektorski moment inercije - warping (net),          Iw = ", "{:.2f}".format(hp1.get_Iw(tc1)), "cm6")
    print(" Polarni moment inercije (net),                      Ip = ", "{:.2f}".format(hp1.get_Ip(tc1)), "cm4")
    print(" St Venant-ov moment inercije (net),                 It = ", "{:.2f}".format(hp1.get_It(tc1)), "cm4")


def SvojstvaHatProfila(h, t, bf, fi, bp, tp):
    #   Provjera karakteristika proizvoljno unesenog Hat profila

    # Eksplicitno zadavanje karakteristika Hat profila
    AH36 = MaterialProperty(1, 210000, 0.3, 7850, 355, "AH36")
    hat1 = HatBeamProperty(1, h, t, bf, fi, AH36)
    tc1 = CorrosionAddition(1, 2)

    # Odabrano sunosivo oplocenje za provjeru

    print(" Visina profila (gross),                             h = ", hat1.h, "mm")
    print(" Debljina struka i prirubnice (gross),               t = ", hat1.t, "mm")
    print(" Sirina prirubnice (gross),                         bf = ", hat1.bf, "mm")
    print(" Kut nagiba struka profila,                         fi = ", hat1.fi, "°")
    print(" Sirina sunosivog oplocenja (gross),                bp = ", bp, "mm")
    print(" Debljina sunosivog oplocenja (gross),              tp = ", tp, "mm")
    print(" Korozijski dodatak,                                tc = ", tc1.tc, "mm")
    print(" Granica razvlacenja,                              Reh = ", hat1.mat.Reh, "N/mm2")
    print(" Razmak strukova Hat profila (net),                 S1 = ", "{:.2f}".format(hat1.getS1_Hat()), "mm")
    print(" Smicna povrsina (net),                           A_sh = ", "{:.2f}".format(hat1.getShArea_Hat(tc1)), "cm2")
    print(" Povrsina samog Hat profila (gross),             A_Hat = ", "{:.2f}".format(hat1.getArea), "cm2")
    print(" Povrsina Hat profila (net) sa oplatom,              A = ", "{:.2f}".format(hat1.getArea_I(bp, tp, tc1)), "cm2")
    print(" Težište Hat profila od simetrale prirubnice,     z_NA = ", "{:.2f}".format(hat1.get_z_na_Hat(bp, tp, tc1)), "mm")
    print(" Moment inercije (net),                             Iy = ", "{:.2f}".format(hat1.get_Iy_I(bp, tp, tc1)), "cm4")
    print(" Moment otpora (net),                             Wmin = ", "{:.2f}".format(hat1.getWmin(bp, tp, tc1)), "cm3")
    print(" Sektorski moment inercije - warping (net),         Iw = ", "{:.2f}".format(hat1.get_Iw(tc1)), "cm6")
    print(" Polarni moment inercije (net),                     Ip = ", "{:.2f}".format(hat1.get_Ip(tc1)), "cm4")
    print(" St Venant-ov moment inercije (net),                It = ", "{:.2f}".format(hat1.get_It(tc1)), "cm4")


def SvojstvaFBProfila(hw, tw, bp, tp):
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    beam = FBBeamProperty(1, hw, tw, ST24)
    tc1 = CorrosionAddition(1, 2)
    print(" Visina struka (gross),                          hw = ", beam.hw, "mm")
    print(" Debljina struka (gross),                        tw = ", beam.tw, "mm")
    print(" Sirina prirubnice (gross),                      bf = ", beam.bf, "mm")
    print(" Debljina prirubnice (gross),                    tf = ", beam.tf, "mm")
    print(" Sirina sunosivog oplocenja,                     bp = ", bp, "mm")
    print(" Debljina sunosivog oplocenja (gross),           tp = ", tp, "mm")
    print(" Korozijski dodatak,                             tc = ", tc1.tc, "mm")
    print(" Granica razvlacenja,                           Reh = ", beam.mat.Reh, "N/mm2")
    print(" Smicna povrsina (net),                        A_sh = ", beam.getShArea_T(tp, tc1), "cm2")
    print(" Povrsina samog T profila (gross),              A_T = ", "{:.2f}".format(beam.getArea), "cm2")
    print(" Povrsina nosaca sa sunosivom sirinom (net),      A = ", "{:.2f}".format(beam.getArea_I(bp, tp, tc1)), "cm2")
    print(" Moment inercije (net),                          Iy = ", "{:.2f}".format(beam.get_Iy_I(bp, tp, tc1)), "cm4")
    print(" Moment otpora (net),                          Wmin = ", "{:.2f}".format(beam.getWmin(bp, tp, tc1)), "cm3")
    print(" Sektorski moment inercije - warping (net),      Iw = ", "{:.2f}".format(beam.get_Iw(tc1)), "cm6")
    print(" Polarni moment inercije (net),                  Ip = ", "{:.2f}".format(beam.get_Ip(tc1)), "cm4")
    print(" St Venant-ov moment inercije (net),             It = ", "{:.2f}".format(beam.get_It(tc1)), "cm4")


def ListaMaterijala(grillage):
    # Spremljena lista materijala
    for i in grillage.material_props().keys():
        mat = grillage.material_props()[i]
        print(" Karakteristike unesenog materijala sa ID =", mat.id)
        print(" Naziv: ", mat.name)
        print(" Modul elasticnosti,     E = ", mat.E, "N/mm2")
        print(" Poissonov koeficijet,   v = ", mat.v)
        print(" Gustoca materijala,    ro = ", mat.ro, "kg/m3")
        print(" Granica razvlacenja,  Reh = ", mat.Reh, "N/mm2", "\n")


def BrojaUnesenihNosaca(grillage):
    print("ID uzduznih nosaca", grillage.longitudinal_members().keys())
    print("ID poprecnih nosaca", grillage.transverse_members().keys())

    print("Broj uzduznih nosaca u longitudinal_members:", len(grillage.longitudinal_members()))
    print("Broj poprecnih nosaca u transverse_members:", len(grillage.transverse_members()))


def JakiNosaci(grillage):
    #   Jaki uzduzni nosaci - smjer, relativna koordinata, koordinate cvorova
    print(" Atributi jakih uzduznih nosaca:")

    # nosac = grillage.longitudinal_members()[1]
    # nosac.flange_direction = FlangeDirection.OUTWARD

    for member in grillage.longitudinal_members().values():
        print("     ID:", member.id, ",", member.direction,
              ", rel. koord. =", "{:.3f}".format(member.rel_dist),
              ", koord u [mm]: node1:", member.end_nodes[0],
              ", node2:", member.end_nodes[1], ", smjer prirubnice:", member.flange_direction)

    #   Jaki poprecni nosaci - smjer, relativna koordinata, koordinate cvorova
    print("\n Atributi jakih poprecnih nosaca:")
    for member in grillage.transverse_members().values():
        print("     ID:", member.id, ",", member.direction,
              ", rel. koord. =", "{:.3f}".format(member.rel_dist),
              ", koordinate u [mm]: node1:", member.end_nodes[0],
              ", node2:", member.end_nodes[1], ", smjer prirubnice:", member.flange_direction)


def CvoroviSegmenata(grillage):
    #   Segmenti svih jakih uzduznih nosaca - koordinate cvorova, duljina segmenta, id nosaca
    print("Segmenti jakih uzduznih nosaca:")
    for i in grillage.longitudinal_members():
        for i_segmenta in range(0, grillage.N_transverse - 1):
            curr_segment = grillage.longitudinal_members()[i].segments[i_segmenta]
            print(" nosac br.", i, ", ID segmenta:", curr_segment.id,
                  ", koordinate node1:", grillage.get_segment_nodes(curr_segment)[0],
                  " , node2:", grillage.get_segment_nodes(curr_segment)[1],
                  ", duljina segmenta:", "{:.3f}".format(grillage.get_segment_length(curr_segment)), "m",
                  ", ID PSM:", curr_segment.primary_supp_mem.id,
                  ", cross1:", curr_segment.cross_member1.id,
                  ", cross2:", curr_segment.cross_member2.id)

    #   Segmenti svih jakih poprecnih nosaca - koordinate cvorova, duljina segmenta,  id nosaca
    print("\n", "Segmenti jakih poprecnih nosaca:")
    for i in grillage.transverse_members():
        for i_segmenta in range(0, grillage.N_longitudinal - 1):
            curr_segment = grillage.transverse_members()[i].segments[i_segmenta]
            print(" nosac br.", i, ", ID segmenta:", curr_segment.id,
                  ", koordinate node1:", grillage.get_segment_nodes(curr_segment)[0],
                  ", node2:", grillage.get_segment_nodes(curr_segment)[1],
                  ", duljina segmenta:", "{:.3f}".format(grillage.get_segment_length(curr_segment)), "m",
                  ", ID PSM:", curr_segment.primary_supp_mem.id,
                  ", cross1:", curr_segment.cross_member1.id,
                  ", cross2:", curr_segment.cross_member2.id)


def DuljineSegmenata(grillage):
    #   Duljina svih segmenata
    print("\n Uzduzni segmenti:")
    for i in grillage.longitudinal_members():
        for i_segmenta in range(0, grillage.N_transverse - 1):
            segment = grillage.longitudinal_members()[i].segments[i_segmenta]
            print("     Segment jakog uzduznog nosaca br.", i, ", ID segmenta:", segment.id,
                  ", Segment length =", "{:.2f}".format(Segment.segment_len(segment)), "m")

    print("\n Poprecni segmenti:")
    for i in grillage.transverse_members():
        for i_segmenta in range(0, grillage.N_longitudinal - 1):
            segment = grillage.transverse_members()[i].segments[i_segmenta]
            print("     Segment jakog poprecnog nosaca br.", i, ", ID segmenta:", segment.id,
                  ", Segment length =", "{:.2f}".format(Segment.segment_len(segment)), "m")


def SvojstvaSegmenata(grillage):
    #   Moment otpora svih segmenata uzduznih nosaca - id nosaca, id segmenta, Wmin, sunosivo oplocenje
    for i in grillage.longitudinal_members():
        for i_segmenta in range(0, grillage.N_transverse - 1):
            curr_segment = grillage.longitudinal_members()[i].segments[i_segmenta]
            print("Jaki uzduzni nosac", i, ", ID segmenta:", curr_segment.id,
                  " ,BeamProperty ID:", curr_segment.beam_prop.id,
                  ", Wmin =", "{:.2f}".format(curr_segment.Wmin), "cm3",
                  ", sunosivo oplocenje: bp =", "{:.1f}".format(Segment.get_attplate(curr_segment)[0]), "mm",
                  ", tp =", Segment.get_attplate(curr_segment)[1], "mm",
                  ", tip profila: ", curr_segment.beam_prop.beam_type,
                  ", tf =", curr_segment.beam_prop.tf,
                  ", tf_net =", curr_segment.beam_prop.tf_net(grillage.corrosion_addition()[1]))
    print("\n")
    for i in grillage.transverse_members():
        for i_segmenta in range(0, grillage.N_longitudinal - 1):
            curr_segment = grillage.transverse_members()[i].segments[i_segmenta]
            print("Jaki poprecni nosac", i, ", ID segmenta:", curr_segment.id,
                  ", BeamProperty ID:", curr_segment.beam_prop.id,
                  ", Wmin =", "{:.2f}".format(curr_segment.Wmin), "cm3",
                  ", sunosivo oplocenje: bp =", "{:.1f}".format(Segment.get_attplate(curr_segment)[0]), "mm",
                  ", tp =", Segment.get_attplate(curr_segment)[1], "mm",
                  ", tip profila: ", curr_segment.beam_prop.beam_type)


def PoljaOplate(grillage):
    #   Sva polja oplate - id, tp, Reh, ID segmenata, nacin zadavanja layouta, vrijednost layouta, orijentacija ukrepa
    for plate_id in range(1, (grillage.N_longitudinal - 1) * (grillage.N_transverse - 1) + 1):
        plate = grillage.plating()[plate_id]
        print(" Polje oplate:", plate.id,
              ",  tp =", plate.plate_prop.tp, "mm",
              ", tp_net =", PlateProperty.tp_net(grillage.corrosion_addition()[1], plate.plate_prop.tp),
              ",  Reh =", plate.plate_prop.plate_mat.Reh,
              # ",  segm.", grillage.plating()[plate_id].long_seg1.id, ",",
              # grillage.plating()[plate_id].long_seg2.id,
              # ",  popr. segm.", grillage.plating()[plate_id].trans_seg1.id, ",",
              # grillage.plating()[plate_id].trans_seg2.id,
              ", def type: ", plate.stiff_layout.definition_type,
              ", iznos: ", plate.stiff_layout.definition_value,
              " ,", plate.stiff_dir)


def DimPoljaOplate(grillage):
    #   Dimenzije svih polja oplate
    for plate_id in range(1, (grillage.N_longitudinal - 1) * (grillage.N_transverse - 1) + 1):
        oplata = grillage.plating()[plate_id]
        print(" Polje oplate:", oplata.id, ", uzduzna dimenzija =", "{:.3f}".format(Plate.plate_longitudinal_dim(oplata)), "m",
              ", poprecna dimenzija =", "{:.3f}".format(Plate.plate_transverse_dim(oplata)), "m")


def SimetricnaPoljaOplate(grillage):
    #   Test simetricnih zona oplate
    grillage.assign_symmetric_plating()

    zona_oplate_ukupno = (grillage.N_longitudinal - 1) * (grillage.N_transverse - 1)
    for plate_id in range(1, zona_oplate_ukupno + 1):
        oplata = grillage.plating()[plate_id]
        if len(oplata.symmetric_plate_zones) == 1:
            print("Zona oplate ID:", oplata.id, ", upisano simetricnih:", len(oplata.symmetric_plate_zones),
                  ",     prvi:", oplata.symmetric_plate_zones[0].id)

        if len(oplata.symmetric_plate_zones) == 2:
            print("Zona oplate ID:", oplata.id, ", upisano simetricnih:", len(oplata.symmetric_plate_zones),
                  ",     prvi:", oplata.symmetric_plate_zones[0].id,
                  ",     drugi:", oplata.symmetric_plate_zones[1].id)

        if len(oplata.symmetric_plate_zones) == 3:
            print("Zona oplate ID:", oplata.id, ", upisano simetricnih:", len(oplata.symmetric_plate_zones),
                  ",     prvi:", oplata.symmetric_plate_zones[0].id,
                  ",     drugi:", oplata.symmetric_plate_zones[1].id,
                  ",     treci:", oplata.symmetric_plate_zones[2].id)


def KoordinateUkrepa(grillage):
    #   Koordinate svih ukrepa na svim poljima oplate
    for oplata in grillage.plating().values():
        print("\n", "Polje oplate: ", oplata.id)
        for i_stiff in range(1, int(oplata.get_stiffener_number()) + 1):
            cvor1 = str(oplata.get_stiff_coords(i_stiff)[0])
            cvor2 = str(oplata.get_stiff_coords(i_stiff)[1])
            string1 = "     Ukrepa br." + str(i_stiff) + ", koordinate 1. cvora:" + cvor1
            string2 = ",         koordinate 2. cvora:" + cvor2
            print(string1, string2)


def SimetricniNosaci(grillage):
    # Provjera simetricnih nosaca
    grillage.assign_symmetric_members()

    print("Uzduzna simetrija:")
    for i_uzduzni in grillage.longitudinal_members().keys():
        print("Uzduzni nosac br.", i_uzduzni, ", postoji simetricni:", grillage.longitudinal_members()[i_uzduzni].has_symmetric_memb)
        if grillage.longitudinal_members()[i_uzduzni].symmetric_member is not None:
            print(" Simetricni mu je nosac br.", grillage.longitudinal_members()[i_uzduzni].symmetric_member.id, "\n")

    print("Poprecna simetrija:")
    for i_poprecni in grillage.transverse_members().keys():
        print("Poprecni nosac br.", i_poprecni, ", postoji simetricni:", grillage.transverse_members()[i_poprecni].has_symmetric_memb)
        if grillage.transverse_members()[i_poprecni].symmetric_member is not None:
            print(" Simetricni mu je nosac br.", grillage.transverse_members()[i_poprecni].symmetric_member.id, "\n")


def TrazilicaPripadnogPoljaOplate(grillage, smjer_segmenta: BeamDirection, id_nosaca, id_segmenta):
    #   Trazilica polja oplate kojemu pripada odabrani segment

    if smjer_segmenta == BeamDirection.LONGITUDINAL:
        segment_uzd = grillage.longitudinal_members()[id_nosaca].segments[id_segmenta]
        for i_plate in grillage.plating().keys():
            plate = grillage.plating()[i_plate]
            test = Plate.test_plate_segment(plate, segment_uzd)
            if test:
                print("Jaki nosac broj:", id_nosaca, "segment ID:", id_segmenta, "definira polje oplate broj", i_plate)

    elif smjer_segmenta == BeamDirection.TRANSVERSE:
        segment_pop = grillage.transverse_members()[id_nosaca].segments[id_segmenta]
        for i_plate in grillage.plating().keys():
            plate = grillage.plating()[i_plate]
            test = Plate.test_plate_segment(plate, segment_pop)
            if test:
                print("Jaki nosac broj:", id_nosaca, "segment ID:", id_segmenta, "definira polje oplate broj", i_plate)


def MasaPoklopca(grillage):
    print("Ukupna masa poklopca:", "{:.2f}".format(grillage.getGrillageMass()), "kg")


def PrimarySuppMemSpacingCheck(spacing, span):
    # Spacing of primary supporting members check, according to:
    #   IACS Common Structural Rules, July 2012, Chapter 9, Section 5, 2. Arrangements, 2.2.3
    """
    :param spacing: Spacing of primary supporting members parallel to the direction of ordinary stiffeners.
    :param span: Span of primary supporting member parallel to the direction of ordinary stiffeners.
    :return: True if spacing of primary supporting members is in compliance with Chapter 9, Section 5, 2. Arrangements, 2.2.3.
    """
    maximum_spacing = 1/3 * np.round(span, 4)
    girder_spacing = np.round(spacing, 4)
    if girder_spacing <= maximum_spacing:
        return True
    elif girder_spacing > maximum_spacing:
        return False


def ProvjeraSpacinga(grillage):
    for plate in grillage.plating().keys():
        oplata = grillage.plating()[plate]
        stiffener_direction = oplata.stiff_dir
        if stiffener_direction == BeamDirection.LONGITUDINAL:
            girder_spacing = Plate.plate_longitudinal_dim(oplata)
            girder_span = grillage.L_overall

            segment1 = oplata.long_seg1
            psm1 = segment1.primary_supp_mem
            print("Jaki nosac ID:", psm1.id, ", segment ID:", segment1.id, ", smjer ukrepa:", stiffener_direction,
                  ", span =", girder_span, "m,   Provjera prema 2.2.3:",
                  PrimarySuppMemSpacingCheck(girder_spacing, girder_span))

        elif stiffener_direction == BeamDirection.TRANSVERSE:
            girder_spacing = Plate.plate_transverse_dim(oplata)
            girder_span = grillage.B_overall
            print("Provjera prema 2.2.3 za sva polja oplate:", PrimarySuppMemSpacingCheck(girder_spacing, girder_span))


def AsocijacijaOplateIjakihNosaca(grillage, nosac_id):
    for plate_id in grillage.plating().keys():
        oplata = grillage.plating()[plate_id]
        jaki_nosac = grillage.longitudinal_members()[nosac_id]
        test_psm = Plate.test_plate_psm(oplata, jaki_nosac)
        print(" Polje oplate ID:", plate_id, ", test:", test_psm)


def AsocijacijaOplateIzmeduUzduznihNosaca(grillage, nosac1_id, nosac2_id):
    for plate_id in grillage.plating().keys():
        oplata = grillage.plating()[plate_id]
        jaki_nosac1 = grillage.longitudinal_members()[nosac1_id]
        jaki_nosac2 = grillage.longitudinal_members()[nosac2_id]
        test_psm = Plate.test_plate_between_psm(oplata, jaki_nosac1, jaki_nosac2)
        print(" Polje oplate ID:", plate_id, ", test:", test_psm)


def DimenzijeElementarnogPanelaOplate(grillage):
    for plate in grillage.plating().values():
        for panel in plate.elementary_plate_panels.values():
            ss, ls = panel.get_elementary_plate_dimensions()
            print(" Polje oplate ID:", plate.id, ", elementarni panel:", panel.id, ", ls=:", ls, "m, ss=", ss, "m",
                  ", edge1:", panel.edge_stiffener_1, ", edge2:", panel.edge_stiffener_2, ", edge3:", panel.edge_stiffener_3,
                  ", edge4:", panel.edge_stiffener_4)
        print("\n")


def ZapisCvorovaPolja(grillage):
    plate_id = 1
    oplata = grillage.plating()[plate_id]
    dim_x = 600                                                 # Dimenzija konacnog elementa (razmak cvorova) po x, [mm]
    dim_y = 625                                                 # Dimenzija konacnog elementa (razmak cvorova) po y, [mm]
    n_x = int((Plate.plate_longitudinal_dim(oplata) * 1000 / dim_x) + 1)    # Broj cvorova po uzduznoj dimenziji
    n_y = int((Plate.plate_transverse_dim(oplata) * 1000 / dim_y) + 1)      # Broj cvorova po poprecnoj dimenziji

    # Prva zona oplate - pocetna
    polje1 = np.zeros((n_y, n_x))
    node_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            polje1[stupac, redak] = node_id
            node_id += 1

    # Druga zona oplate - koristi dio cvorova pocetne
    polje2 = np.zeros((n_y, n_x))
    for stupac in range(0, n_y):
        for redak in range(1, n_x):
            polje2[stupac, 0] = polje1[stupac, n_x-1]
            polje2[stupac, redak] = node_id
            node_id += 1

    print(polje1)
    print(polje2)


def ZapisCvorovaOplate(grillage):
    L = grillage.L_overall * 1000
    B = grillage.B_overall * 1000

    dim_x = 600                                                 # Dimenzija konacnog elementa (razmak cvorova) po x, [mm]
    dim_y = 625                                                 # Dimenzija konacnog elementa (razmak cvorova) po y, [mm]
    n_x = int(L / dim_x)    # Broj cvorova po uzduznoj dimenziji
    n_y = int(B / dim_y)      # Broj cvorova po poprecnoj dimenziji

    polje1 = np.zeros((n_y, n_x))

    # Cvorovi cijele oplate
    node_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            polje1[stupac, redak] = node_id
            node_id += 1
    print(polje1)


def TestPlateCommonSegment(grillage, plate_id1, plate_id2):
    plate1 = grillage.plating()[plate_id1]
    plate2 = grillage.plating()[plate_id2]
    segment = Grillage.plate_common_segment(plate1, plate2)
    if segment is None:
        print("Ne postoji zajednicki segment polja oplate", plate_id1, "i", plate_id2)
    else:
        print("Zajednicki segment polja", plate_id1, "i", plate_id2, "ima ID:", segment.id,
              ", koordinate cvorova:", Segment.get_segment_node1(segment), Segment.get_segment_node2(segment))


def TestALlPlateCommonSegments(grillage):
    for plate_id_n in grillage.plating().keys():
        for plate_id_m in grillage.plating().keys():
            if plate_id_m != plate_id_n:
                oplata_n = grillage.plating()[plate_id_n]
                oplata_m = grillage.plating()[plate_id_m]
                comm_segment = Grillage.plate_common_segment(oplata_n, oplata_m)
                if comm_segment is not None:
                    print("Polje:", oplata_n.id, "i polje", oplata_m.id, "imaju zajednicki segment ID:", comm_segment.id,
                          ", nosaca ID:", comm_segment.primary_supp_mem.id, "orijentacije:", comm_segment.primary_supp_mem.direction)


def TestStiffenerSpacingNumber(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    spacing = Plate.get_stiffener_spacing(plate)
    number = Plate.get_stiffener_number(plate)
    offset = Plate.get_equal_stiffener_offset(plate)
    path = Plate.get_path_length(plate)
    print("Spacing:", spacing, "Number:", number, "Offset", offset, "path_length", path)


def TestGetAdjacentPlates(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    adjacent_plates_list = grillage.get_adjacent_plates(plate)
    if adjacent_plates_list is not None:
        for i in range(0, len(adjacent_plates_list)):
            print(adjacent_plates_list[i].id)


def TestGetLongitudinalAdjacentPlates(grillage, plate_id):
    # Susjedna polja u uzdužnom smjeru imaju zajednički poprečni segment
    plate = grillage.plating()[plate_id]
    adjacent_plates_list = grillage.get_long_adjacent_plates(plate)
    if adjacent_plates_list is not None:
        for i in range(0, len(adjacent_plates_list)):
            print("Polju oplate", plate_id, "je u uzdužnom smjeru susjedno polje oplate", adjacent_plates_list[i].id)


def TestGetTransverseAdjacentPlates(grillage, plate_id):
    # Susjedna polja u poprečnom smjeru imaju zajednički uzdužni segment
    plate = grillage.plating()[plate_id]
    adjacent_plates_list = grillage.get_tran_adjacent_plates(plate)
    if adjacent_plates_list is not None:
        for i in range(0, len(adjacent_plates_list)):
            print("Polju oplate", plate_id, "je u poprečnom smjeru susjedno polje oplate", adjacent_plates_list[i].id)


def Test_plating_zones_between_psm(grillage):
    for uzduzni in range(1, grillage.N_longitudinal):
        nosac1 = grillage.longitudinal_members()[uzduzni]
        nosac2 = grillage.longitudinal_members()[uzduzni + 1]
        lista = grillage.plating_zones_between_psm(nosac1, nosac2)
        for i in range(0, len(lista)):
            print("Između uzdužnih nosača", nosac1.id, "i", nosac2.id, "je polje oplate", lista[i].id)

    for poprecni in range(1, grillage.N_transverse):
        nosac1 = grillage.transverse_members()[poprecni]
        nosac2 = grillage.transverse_members()[poprecni + 1]
        lista = grillage.plating_zones_between_psm(nosac1, nosac2)
        for i in range(0, len(lista)):
            print("Između poprečnih nosača", nosac1.id, "i", nosac2.id, "je polje oplate", lista[i].id)


def Test_assign_symmetric_segments(grillage):
    # Provjera simetricnih nosaca
    grillage.assign_symmetric_segments()

    print("Uzduzni nosaci:")
    for i_uzduzni in grillage.longitudinal_members().keys():
        for i_segment in grillage.longitudinal_members()[i_uzduzni].segments:
            if i_segment.symmetric_segment is not None:
                print("Nosac br.", i_uzduzni, " segmentu", i_segment.id, "je simetricni segment", i_segment.symmetric_segment.id, "\n")

    print("Poprecni nosaci:")
    for i_poprecni in grillage.transverse_members().keys():
        for i_segment in grillage.transverse_members()[i_poprecni].segments:
            if i_segment.symmetric_segment is not None:
                print("Nosac br.", i_poprecni, " segmentu", i_segment.id, "je simetricni segment", i_segment.symmetric_segment.id, "\n")


def Test_set_long_symm_segment_beam_property(grillage, prim_supp_member_id: int, segment_id: int):
    grillage.assign_symmetric_members()
    grillage.assign_symmetric_segments()

    beam_property = grillage.beam_props()[5]
    Grillage.set_long_symm_segment_beam_property(grillage, prim_supp_member_id, segment_id, beam_property)

    for i in grillage.longitudinal_members():
        for i_segmenta in range(0, grillage.N_transverse - 1):
            curr_segment = grillage.longitudinal_members()[i].segments[i_segmenta]
            beam = curr_segment.beam_prop
            print("Jaki uzduzni nosac", i, ", ID segmenta:", curr_segment.id,
                  " ,BeamProperty ID:", curr_segment.beam_prop.id,
                  ", Wmin =", "{:.2f}".format(curr_segment.Wmin), "cm3",
                  ", hw =", beam.hw, "mm",
                  ", tw =", beam.tw, "mm",
                  ", bf =", beam.bf, "mm",
                  ", tf =", beam.tf, "mm")


def Test_set_long_member_beam_property(grillage, prim_supp_member_id: int):
    beam_property = grillage.beam_props()[5]    # Novi beam property
    Grillage.set_long_member_beam_property(grillage, prim_supp_member_id, beam_property)

    for i in grillage.longitudinal_members():
        for i_segmenta in range(0, grillage.N_transverse - 1):
            curr_segment = grillage.longitudinal_members()[i].segments[i_segmenta]
            beam = curr_segment.beam_prop
            print("Jaki uzduzni nosac", i, ", ID segmenta:", curr_segment.id,
                  " ,BeamProperty ID:", curr_segment.beam_prop.id,
                  ", Wmin =", "{:.2f}".format(curr_segment.Wmin), "cm3",
                  ", hw =", beam.hw, "mm",
                  ", tw =", beam.tw, "mm",
                  ", bf =", beam.bf, "mm",
                  ", tf =", beam.tf, "mm")


def Test_Wmin_Iy_ukrepa(grillage):
    for plate in grillage.plating().values():
        print("Zona oplate", plate.id, ", Karakteristike ukrepa: Wmin =", "{:.2f}".format(plate.Wmin()),
              "cm3, Iy =", "{:.2f}".format(plate.Iy()), "cm4")


def Test_segments_between_psm(grillage, nosac1_id, nosac2_id, direction: BeamDirection):
    if direction == BeamDirection.LONGITUDINAL:
        nosac1 = grillage.longitudinal_members()[nosac1_id]
        nosac2 = grillage.longitudinal_members()[nosac2_id]
    else:
        nosac1 = grillage.transverse_members()[nosac1_id]
        nosac2 = grillage.transverse_members()[nosac2_id]

    segment_list = grillage.segments_between_psm(nosac1, nosac2)
    print("Upisano segmenata:", len(segment_list))
    print(segment_list)


def Test_midpoint():
    # Koordinate središta između proizvoljno odabranih tocaka

    # Točke:
    # node1 = Node(1, 0, 6160, 1200)
    # node2 = Node(2, 0, 12020, 1200)

    node1 = ModelNode(1, 0, 12020, 1200)
    node2 = ModelNode(2, 6430, 6160, 1200)

    # Vracanje koordinata kroz metodu coords
    n1 = node1.coords
    n2 = node2.coords

    midpoint = Grillage.get_midpoint(n1, n2)

    print("Unesene točke:", n1, ",", n2)
    print("Koordinate središta", midpoint)
    print(type(midpoint))


def Test_get_elementary_plate(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    for panel in plate.elementary_plate_panels.values():
        ss, ls = panel.get_elementary_plate_dimensions()
        n_sub = panel.sub_panel_number
        if panel.stiffener_1_id is not None:
            stiff_1_coords = plate.get_stiff_coords(panel.stiffener_1_id)
        else:
            stiff_1_coords = None
        if panel.stiffener_2_id is not None:
            stiff_2_coords = plate.get_stiff_coords(panel.stiffener_2_id)
        else:
            stiff_2_coords = None

        print(" Polje oplate ID:", plate.id, ", elementarni panel:", panel.id)
        print("     Dimenzije sub panela za proracun podobnosti:   ls=:", "{:.3f}".format(ls), "m, ss=", "{:.3f}".format(ss), "m")
        print("     Broj interkostalnih ukrepa: ", panel.intercostal_stiffener_num,  ", broj sub panela:", n_sub)
        print("     ID prve ukrepe:", panel.stiffener_1_id,  ", ID druge ukrepe:", panel.stiffener_2_id)
        print("     Koordinate: prve ukrepe:", stiff_1_coords,  ", druge ukrepe:", stiff_2_coords)

        for sub_panel in range(1, n_sub + 1):
            edge_stiff = panel.get_edge_beam_types(sub_panel)
            print("     Rubni nosaci sub panela", sub_panel, ":   ", edge_stiff)
        print("\n")


def Test_intercostal_coords(grillage, plate_id, elementary_panel_id, intercostal_n):
    plate = grillage.plating()[plate_id]
    elementary_panel = plate.elementary_plate_panels[elementary_panel_id]
    print("Koordinata čvora interkostalne ukrepe na zoni oplate", plate_id, ", elementarni panel broj", elementary_panel_id, ", ukrepa", intercostal_n)
    elementary_panel.get_intercostal_coords(intercostal_n)


def Test_all_intercostal_coords(grillage):
    for plate in grillage.plating().values():
        print("Zona oplate", plate.id)
        for elementary_plate in plate.elementary_plate_panels.values():
            num_of_intercostals = elementary_plate.intercostal_stiffener_num
            if num_of_intercostals == 0:
                print("  nema interkostalnih ukrepa")
                break

            print("  Elementarni panel broj", elementary_plate.id)
            for intercostal in range(1, num_of_intercostals + 1):
                node1, node2 = elementary_plate.get_intercostal_coords(intercostal)
                print("     Interkostalna ukrepa broj", intercostal, ", koordinate:", node1, ",", node2)


def PlotGrillageTopology(grillage):
    import matplotlib.pyplot as plt

    # Plot svih segmenata jakih uzduznih nosaca
    for i_longitudinal in grillage.longitudinal_members():
        for i_segment in range(0, grillage.N_transverse - 1):
            curr_segment = grillage.longitudinal_members()[i_longitudinal].segments[i_segment]
            nodes = grillage.get_segment_nodes(curr_segment)
            x1 = nodes[0][0]
            x2 = nodes[1][0]
            y1 = nodes[0][1]
            y2 = nodes[1][1]
            plt.plot([x1, x2], [y1, y2], "g", linewidth=2)  # "g" daje zelenu boju za jake uzduzne nosace

    # Plot svih segmenata jakih poprecnih nosaca
    for i_transverse in grillage.transverse_members():
        for i_segment in range(0, grillage.N_longitudinal - 1):
            curr_segment = grillage.transverse_members()[i_transverse].segments[i_segment]
            nodes = grillage.get_segment_nodes(curr_segment)
            x1 = nodes[0][0]
            x2 = nodes[1][0]
            y1 = nodes[0][1]
            y2 = nodes[1][1]
            plt.plot([x1, x2], [y1, y2], "r", linewidth=2)  # "r" daje crvenu boju za jake poprecne nosace
            #                                                   c=np.random.rand(3,) daje nasumicnu boju

    # Plot svih ukrepa
    for id_oplate in grillage.plating().keys():
        oplata = grillage.plating()[id_oplate]
        num_of_stiff = int(Plate.get_stiffener_number(oplata))
        for i_stiff in range(1, num_of_stiff + 1):
            nodes = Plate.get_stiff_coords(oplata, i_stiff)
            x1 = nodes[0][0]
            x2 = nodes[1][0]
            y1 = nodes[0][1]
            y2 = nodes[1][1]
            plt.plot([x1, x2], [y1, y2], "--k", linewidth=1.5)  # "--k" daje crnu isprekidanu liniju

    # Plot svih interkostalnih ukrepa
    for plate in grillage.plating().values():
        for elementary_plate in plate.elementary_plate_panels.values():
            num_of_intercostals = elementary_plate.intercostal_stiffener_num
            for intercostal in range(1, num_of_intercostals + 1):
                nodes = elementary_plate.get_intercostal_coords(intercostal)
                x1 = nodes[0][0]
                x2 = nodes[1][0]
                y1 = nodes[0][1]
                y2 = nodes[1][1]
                plt.plot([x1, x2], [y1, y2], "--k", linewidth=0.8)  # "--k" daje crnu isprekidanu liniju

    plt.ylim(-1000, grillage.B_overall * 1000 + 1000)  # Granice plota po y
    plt.xlim(-1000, grillage.L_overall * 1000 + 1000)  # Granice plota po x
    plt.gca().set_aspect('equal')  # Fix za aspect ratio
    plt.show()


#   ******************   1. Dio: Lista opcih provjera metoda   ******************

# KoordinateCvorova()
# SvojstvaTprofila(1001, 11, 151, 21, 500, 9)
# SvojstvaHPprofila(200, 11, 404, 9)
# SvojstvaHatProfila(150, 10, 200, 50, 600, 10)
# SvojstvaHatProfila(220, 10, 220, 80, 600, 10)
# SvojstvaFBProfila(1000, 12, 0, 0)

#   ******************   2. Dio: Lista provjera za odabranu topologiju   ******************
# ListaMaterijala(hc_variant)
# BrojaUnesenihNosaca(hc_variant)
# JakiNosaci(hc_variant)
# CvoroviSegmenata(hc_variant)
# DuljineSegmenata(hc_variant)
# SvojstvaSegmenata(hc_variant)
# PoljaOplate(hc_variant)
# DimPoljaOplate(hc_variant)
# SimetricnaPoljaOplate(hc_variant)
# KoordinateUkrepa(hc_variant)
# SimetricniNosaci(hc_variant)
# TrazilicaPripadnogPoljaOplate(hc_variant, BeamDirection.TRANSVERSE, 2, 3)
# MasaPoklopca(hc_variant)
# ProvjeraSpacinga(hc_variant)
# AsocijacijaOplateIjakihNosaca(hc_variant, 1)
# AsocijacijaOplateIzmeduUzduznihNosaca(hc_variant, 1, 2)
# TestPlateCommonSegment(hc_variant, 2, 5)
# TestALlPlateCommonSegments(hc_variant)
# ZapisCvorovaPolja(hc_variant)
# ZapisCvorovaOplate(hc_variant)
# TestStiffenerSpacingNumber(hc_variant, 1)
# TestGetAdjacentPlates(hc_variant, 2)
# TestGetLongitudinalAdjacentPlates(hc_variant, 2)
# TestGetTransverseAdjacentPlates(hc_variant, 6)
# Grillage.hc_variant_check(hc_variant)
# Test_plating_zones_between_psm(hc_variant)
# Test_assign_symmetric_segments(hc_variant)
# Test_set_long_symm_segment_beam_property(hc_variant, 3, 1)
# Test_set_long_member_beam_property(hc_variant, 3)
# Test_Wmin_Iy_ukrepa(hc_variant)
# Test_segments_between_psm(hc_variant, 2, 3, BeamDirection.LONGITUDINAL)
# Test_midpoint()
# DimenzijeElementarnogPanelaOplate(hc_variant)
# Test_get_elementary_plate(hc_variant, 6)
# Test_intercostal_coords(hc_variant, 1, 1, 1)
# Test_all_intercostal_coords(hc_variant)
# PlotGrillageTopology(hc_variant)

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
