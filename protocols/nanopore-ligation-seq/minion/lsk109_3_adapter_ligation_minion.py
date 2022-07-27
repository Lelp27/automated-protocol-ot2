from opentrons import types, protocol_api

# this data is from Ligation Sequencing Kit: DNA Repair and End-Prep 
# modified script from Sakib <sakib.hossain@opentrons.com>
metadata = {
    'protocolName': 'Nanopore LSK109(3) Adapter_Ligation (AMX)',
    'author': 'Seong-Kun Bak <tjdrns27@kribb.re.kr>',
    'apiLevel': '2.11',
    'description': 'with movable PCR tube pause the protocol when the tube should move'
}

def run(protocol: protocol_api.ProtocolContext):
    #============= Deck setting

    ## Module setting
    ### for the movement of Pippete Input all of modules in Device!
    module_magnetic = protocol.load_module('magnetic module gen2', '4')
    module_thermocycler = protocol.load_module("thermocycler Module")

    ## Labware setting
    enzyme_rack = protocol.load_labware("opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap", 1)
    falcon_rack = protocol.load_labware("opentrons_6_tuberack_falcon_50ml_conical", 2)
    dna_plate = protocol.load_labware("bioneer_96_tuberack_200ul", 5)
    mag_plate = module_magnetic.load_labware("biorad_96_tuberack_200ul", 4)
    trash = protocol.loaded_labwares[12]["A1"]

    tiprack_20_sin = protocol.load_labware("opentrons_96_tiprack_20ul", 3)
    tiprack_300_sin = protocol.load_labware("opentrons_96_tiprack_300ul", 6)
    p20_sin = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[tiprack_20_sin])
    p300_sin = protocol.load_instrument("p300_single_gen2", "left", tip_racks=[tiprack_300_sin])

    ## Start Tiprack positions
    p20_sin.starting_tip = tiprack_20_sin.well("A4")
    p300_sin.starting_tip = tiprack_300_sin.well("A1")

    ## Reagents
    ### All of reagents in Eppendorf 1.5ml lobind tube.
    LNB = enzyme_rack["A1"]
    T4_ligase = enzyme_rack["A2"]
    AMX = enzyme_rack["A3"]
    LFB_SFB = enzyme_rack["A4"]
    EB = enzyme_rack["A5"]

    ampure_beads = falcon_rack["A1"]
    

    ## Sample Wells
    start_position = 0 ## start from A1
    sample_number = 2

    dna_well = [dna_plate.wells()[i] for i in range(start_position, start_position + sample_number)]
    final_well = [dna_plate.wells()[i] for i in range(start_position + 16, start_position + 16 + sample_number)]
    # magnetic의 A1이 위치상 자석과 거리가 멀다..
    mag_well = [mag_plate.wells()[i] for i in range(start_position, start_position + sample_number)]
    

    #============= PROTOCOL STEPS
    ## Functions
    def enzyme_transfer(pipette, volume, src, dest, delay_second, asp_rate=None, dis_rate=None, mix_after=False):
        if asp_rate != None:
            pipette.flow_rate.aspirate=asp_rate
        if dis_rate != None:
            pipette.flow_rate.dispense=dis_rate            
        pipette.pick_up_tip()
        pipette.aspirate(volume, src)
        protocol.delay(seconds=delay_second)
        pipette.touch_tip()
        pipette.dispense(volume, dest)
        protocol.delay(seconds=delay_second/2)
        if type(mix_after) == tuple:
            try:
                pipette.flow_rate.aspirate=mix_after[2]
                pipette.flow_rate.dispense=mix_after[2]
            except:
                pass
            pipette.mix(mix_after[0], mix_after[1])
        pipette.blow_out()    
        pipette.touch_tip()
        pipette.drop_tip()

    def remove_supernantant(pipette, volume, src, dest, asp_rate=20, z=0):
        pipette.pick_up_tip()
        pipette.flow_rate.aspirate=asp_rate
        pipette.aspirate(volume, src.bottom().move(types.Point(x=0, y=0, z=z))) #biorad tuberack
        pipette.dispense(volume, dest)
        pipette.drop_tip()


    for dest in dna_well:
        enzyme_transfer(p300_sin, 25, LNB, dest, 1, 25, 25)
    for dest in dna_well:
        enzyme_transfer(p20_sin, 10, T4_ligase, dest, 1, 10, 10)
    for dest in dna_well:
        enzyme_transfer(p20_sin, 5, AMX, dest, 1, 5, 5, mix_after=(3,20))
    # mix 속도가 너무느림

    protocol.pause("1. Check enzyme\
            \n2. Spindown\
            \n3. Resume")

    protocol.delay(minutes=10)

    ## Resuspend AMPure Beads (4)
    p300_sin.flow_rate.aspirate=300
    p300_sin.flow_rate.dispense=300
    p300_sin.pick_up_tip()
    p300_sin.mix(5, 300, ampure_beads.bottom(z=15))
    p300_sin.blow_out()
    p300_sin.move_to(ampure_beads.top())
    protocol.delay(seconds=1)
    p300_sin.touch_tip()
    p300_sin.drop_tip()

    for dest in dna_well:
        enzyme_transfer(p300_sin, 50, ampure_beads, dest, 1.5, 50, 50, mix_after=(3, 30))
    

    ## transfer to tube for Hula Mixer (5)
    protocol.pause("1. Move sample tube to Hula Mixer\
                \n2. Put sample tube in Magnetic well\
                \n3. Resume")

    ## Engage Magnet and Delay for 5 Minutes (6)
    module_magnetic.engage(height_from_base=3.5) # Biorad tube rack = 3
    protocol.delay(minutes=8)

    ## Remove Supernant (7)
    for src in mag_well:
        remove_supernantant(p300_sin, 130, src, trash, asp_rate=50, z=0)

    ## LFB_SFB wash
    # 아직 검증못함.
    module_magnetic.disengage()
    for _ in range(1):
        
        for src in mag_well:
            p300_sin.pick_up_tip()
            p300_sin.flow_rate.aspirate=200
            p300_sin.flow_rate.dispense=150
            p300_sin.aspirate(200, LFB_SFB)
            p300_sin.touch_tip()
            p300_sin.dispense(200, src.top().move(types.Point(2, 0, -15)))
            p300_sin.mix(3, 200, src.top().move(types.Point(2, 0, -15)), rate=1.5)
            p300_sin.move_to(src.top())
            p300_sin.blow_out()
            p300_sin.touch_tip()
            p300_sin.drop_tip()
        ## Remove

    for src in mag_well:
        remove_supernantant(p300_sin, 220, src, trash, asp_rate=100, z=0)
        #remove시 액체가 너무 많아 넘친다(피펫이천천히 내려가야 할 것 같음 or 볼륨을 줄이기)

    ## Pause and Remove Samples for Spin Down (15)
    module_magnetic.engage(height_from_base=2)
    protocol.pause('1. Spindown\
                \n2. replace tube in magnetic_well\
                \n3. Resume')

    ## Remove Residual Ethanol (16)
    protocol.delay(seconds=30)
    for src in mag_well:
        remove_supernantant(p20_sin, 15, src, trash, asp_rate=10, z=-0.5)
    module_magnetic.disengage()
    protocol.delay(seconds=40)

    module_thermocycler.open_lid()

    ## Add 15 of Elution buffer
    for dest in mag_well:
        enzyme_transfer(p20_sin, 15, EB, dest, 1, 15, 15, mix_after=(3, 10))
    protocol.pause("1. If you want get high-molecular weight, move tube to thermocycler \
                  \n2. or Skip this step now")
    module_thermocycler.set_block_temperature(37)
    module_thermocycler.set_lid_temperature(37)
    protocol.delay(minutes=10)

    module_magnetic.engage(height_from_base=2.5)
    protocol.delay(minutes=1)
    protocol.pause("Final step : Elution")
    module_magnetic.disengage()