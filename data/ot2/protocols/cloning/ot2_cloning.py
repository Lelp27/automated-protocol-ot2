# 23-10-20 Simulated by kun.

from opentrons import protocol_api
# from opentrons import simulate
# protocol = simulate.get_protocol_api("2.13")
import pandas as pd
import re
import time
import json
import requests

metadata = {
    "protocolName": "{{PRESENT_TIME}} Cloning (PCR, Assembly, Transformation)",
    "author": "Seong-Kun Bak <sanekun@kaist.ac.kr>",
    "apiLevel": "2.13",
    "robotType": "OT-2",
    "description": "Cloning Protocol in SBL",
}

# [PARAMETERS]
PARAMETERS = {{EXPORT_JSON}}
debug = False

# [Functions]
def discord_message(message, user=PARAMETERS["Parameters"]["Messenger"]):
    if user == "None":
        return None

    # Send message to discord chaneel
    url = "https://discord.com/api/webhooks/1169608986891390996/wLti5ulSOXBtGelwiYO4SJi1jQSJ9tEQ8uC8sVXVlWXih3XWwkoLDr6cJBEm2iPY9b0t"

    headers = {"Content-Type": "application/json"}
    data = {"content": message}
    response = requests.post(url, headers=headers, data=json.dumps(data), verify=False)

def flow_rate(pipette, **kwargs):
    # Change flow rate of pipette

    assert (
        item in ["aspirate", "dispense", "blow_out"] for item in kwargs.keys()
    ), "Error Keywords in Flow Rate."
    for i in kwargs.keys():
        setattr(pipette.flow_rate, i, kwargs[i])

def find_materials_well(material, type: ["DNA", "Enzyme"],
                        PARAMETERS=PARAMETERS, right_well=False):
    # Convert material name to well
    # right_well means well of right side of plate (for abstraction)

    if type == "DNA":
        for plate in PARAMETERS["Plates"].values():
            if plate["type"] != "TF":
                for well in plate["data"].keys():
                    if plate["data"][well] == material:
                        wells = plate["Deck"].wells_by_name()
                        if right_well:
                            pos = list(wells.keys()).index(well)
                            return wells[list(wells)[pos+8]]
                        else:
                            return wells[well]
    if type == "Enzyme":
        return PARAMETERS["Deck"]["Enzyme_position"][material]
    if type == "DW":
        return PARAMETERS["Deck"]["Enzyme_position"][material]


def transfer_materials(key, p20, p300, mix_last=(0,0)):
    # Transfer materials by type of materials

    reaction = PARAMETERS["Reactions"][key]
    volume_dict = PARAMETERS["Reaction_volume"][reaction["type"]]

    # transform dict to dataframe (row iterable)
    df = pd.DataFrame(reaction["data"])

    # Only First Material transfer to all wells at once (Enzyme1 or DW)
    for dw in df["DW"].unique():
        tmp = df[df["DW"] == dw]
        # If Data doens't exist, this step will be skipped
        if len(tmp):
            if pd.isna(dw) or dw == "":
                continue
            src = find_materials_well(dw, "DW")
            vol = float(volume_dict["DW"])
            dest = [find_materials_well(name, "DNA") for name in tmp["Name"].values]

            flow_rate(p300, aspirate=50, dispense=50, blow_out=20)
            p300.distribute(
                vol, src, dest,
                new_tip="once", touch_tip=False,
                disposal_volume=5,
                blow_out=False,
                trash=not debug
            )

    for enzyme_name in df["Enzyme1"].unique():
        tmp = df[df["Enzyme1"] == enzyme_name]
        # If Data doens't exist, this step will be skipped
        if len(tmp):
            if pd.isna(enzyme_name) or enzyme_name == "":
                continue
            src = find_materials_well(enzyme_name, "Enzyme").bottom(z=3)
            vol = float(volume_dict["Enzyme1"])
            dest = [find_materials_well(name, "DNA").bottom(z=1) for name in tmp["Name"].values]

            flow_rate(p300, aspirate=20, dispense=20, blow_out=20)
            p300.distribute(
                vol, src, dest,
                new_tip="once", touch_tip=False,
                mix_before=(2, 50), disposal_volume=5,
                blow_out=True, blowout_location="source well",
                trash=not debug
            )
    # Other Materials
    df.drop(columns=["Enzyme1", "DW"], inplace=True)
    columns = df.columns
    for action in df.values:
        for sample_type, sample_name in zip(columns, action):
            # Empty well will be skipped
            if sample_type == "Name":
                dest = find_materials_well(sample_name, "DNA")
                continue
            if pd.isna(sample_name) or sample_name == "":
                continue

            src = find_materials_well(sample_name, re.sub(r"[0-9]+", "", sample_type))
            vol = float(volume_dict[sample_type])

            # Check DNA or Enzyme
            if re.sub(r"[0-9]+", "", sample_type) == "Enzyme":
                flow_rate(p20, aspirate=1, dispense=1, blow_out=1)
                touch_tip = False
            else:  # DNA & DW
                flow_rate(p20, aspirate=1, dispense=1, blow_out=1)
                touch_tip = False

            p20.transfer(
                vol, src, dest, new_tip="always",
                touch_tip=touch_tip,
                blow_out=False, blowout_location="destination well",
                trash=not debug
            )

        # Mix Product
        if sum(mix_last):
            flow_rate(p20, aspirate=10, dispense=10, blow_out=10)
            p20.pick_up_tip()
            for _ in range(mix_last[0]):
                p20.aspirate(mix_last[1], dest.bottom())
                p20.dispense(mix_last[1], dest.bottom(z=3))
            p20.drop_tip()


def spotting_dispense(pipette, src, dest: list, spotting_volume=4):
    # Dispense liquid at the top of well
    # and then, move down pipette to specific height.

    if not pipette.has_tip:
        pipette.pick_up_tip()

    disposal_vol = 1
    whole_vol = spotting_volume * len(dest) + disposal_vol
    cnt = 0
    while (whole_vol > spotting_volume):
        if pipette.max_volume < whole_vol:
            pipette.aspirate(pipette.max_volume, src)
        else:
            pipette.aspirate(whole_vol, src)

        while (pipette.current_volume > spotting_volume):
            pipette.dispense(spotting_volume, dest[cnt].bottom(z=4.4))
            pipette.move_to(dest[cnt].bottom(z=3))
            cnt += 1
            whole_vol -= spotting_volume

        pipette.blow_out(pipette.trash_container.wells()[0])


def run(protocol: protocol_api.ProtocolContext):
    discord_message(f"Protocol Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Deck Setting
    ## Modules
    tc_mod = protocol.load_module(module_name="thermocyclerModuleV1")
    tc_mod.open_lid()
    Enzyme_deck = protocol.load_labware("opentrons_24_tuberack_nest_1.5ml_screwcap", 1)

    ## Pipette
    p20_tip = protocol.load_labware("opentrons_96_tiprack_20ul", 2)
    p300_tip = protocol.load_labware("opentrons_96_tiprack_300ul", 3)
    p20 = protocol.load_instrument("p20_single_gen2", "left", tip_racks=[p20_tip])
    p300 = protocol.load_instrument("p300_single_gen2", "right", tip_racks=[p300_tip])

    ## Enzymes
    for key in PARAMETERS["Deck"]["Enzyme_position"].keys():
        PARAMETERS["Deck"]["Enzyme_position"][key] = Enzyme_deck.wells_by_name()[
            PARAMETERS["Deck"]["Enzyme_position"][key]
        ]

    ## Plates
    for key in PARAMETERS["Plates"].keys():
        location = PARAMETERS["Deck"]["Deck_position"][key]
        if str(location) == "7":
            PARAMETERS["Plates"][key]["Deck"] = tc_mod.load_labware(
                PARAMETERS["Plates"][key]["labware"])
            continue

        PARAMETERS["Plates"][key]["Deck"] = protocol.load_labware(
            PARAMETERS["Plates"][key]["labware"], location=location)

    # RUN first Reaction
    # PCR
    reaction = list(PARAMETERS["Reactions"].keys())[0]
    final_volume = sum([float(i) for i in PARAMETERS["Reaction_volume"]["PCR"].values()])
    Run_PCR = [i for i in PARAMETERS["Reactions"][reaction]["data"]["Name"].values() if i]

    if Run_PCR:
        discord_message(f"PCR Start")
        transfer_materials(key=reaction, p20=p20, p300=p300, mix_last=(2, 15))
        ## Thermocycler
        discord_message(f"Thermocycler in PCR start RUN take off Enzyme")
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(94)
        tc_mod.set_block_temperature(
            temperature=94,
            hold_time_seconds=30,
            block_max_volume=final_volume
        )
        profile = [
            {"temperature": 94, "hold_time_seconds": 20},
            {"temperature": 55, "hold_time_seconds": 20},
            {"temperature": 68, "hold_time_seconds": int(PARAMETERS["Parameters"]["PCR_extension_time"])},
        ]
        tc_mod.execute_profile(
            steps=profile,
            repetitions=30,
            block_max_volume=final_volume
        )
        tc_mod.set_block_temperature(
            temperature=68,
            hold_time_seconds=60,
            block_max_volume=final_volume
        )
        discord_message(f"PCR will be end 5 minutes later, Take in Enzyme for next step")
        tc_mod.set_block_temperature(
            temperature=8,
            hold_time_minutes=5,
            block_max_volume=final_volume
        )
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(
            temperature = 8,
            hold_time_minutes=2,
            block_max_volume=final_volume
            )
        tc_mod.open_lid()

    # Second Reaction
    reaction = list(PARAMETERS["Reactions"].keys())[1]
    final_volume = sum([float(i) for i in PARAMETERS["Reaction_volume"]["PCR"].values()])
    Run_Assembly = [i for i in PARAMETERS["Reactions"][reaction]["data"]["Name"].values() if i]

    ## If not Empty Run protocol
    if Run_Assembly:
        if PARAMETERS["Parameters"]["Stop_between_reactions"]:
            protocol.pause(f"Place down {list(PARAMETERS['Reactions'].keys())[1]} enzyme")
        discord_message(f"Assembly Start")

        transfer_materials(list(PARAMETERS["Reactions"].keys())[1], p20=p20, p300=p300, mix_last=(2, 10))
        discord_message(f"Thermocycler in Assembly start RUN remove Enzyme")
        ## Thermocycler
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(80)

        ### DpnI
        tc_mod.set_block_temperature(
            temperature=37,
            hold_time_minutes=5,
            block_max_volume=final_volume
        )
        ### denaturation
        tc_mod.set_block_temperature(
            temperature=65,
            hold_time_seconds=20,
            block_max_volume=final_volume
        )
        tc_mod.set_block_temperature(
            temperature=50,
            hold_time_minutes=30,
            block_max_volume=final_volume,
        )

        discord_message(f"Assembly will be end 10 minutes later, Take in CP cell for next step")
        # Ramp rate is 0.1 degree per second
        current_tmp = 50
        while True:
            current_tmp -= 5

            if current_tmp <= 8:
                tc_mod.deactivate_lid()
                tc_mod.set_block_temperature(
                    temperature = 8,
                    hold_time_minutes=2,
                    block_max_volume=final_volume
                    )
                tc_mod.open_lid()
                break

            tc_mod.set_block_temperature(
                temperature=current_tmp,
                hold_time_seconds=45, # ramp_rate is almost 0.1 degree per second
                block_max_volume=final_volume
            )

    # Transformation
    ## IF Tf plate not empty, run protocol
    Run_TF=False
    for key in PARAMETERS["Plates"].keys():
        plate = PARAMETERS["Plates"][key]
        if plate["type"] == "TF":
            if plate["data"]:
                Run_TF=True
                break

    if Run_TF:
        discord_message(f"Transformation Start")
        if PARAMETERS["Parameters"]["Stop_between_reactions"]:
            protocol.pause("Place down TF-associated materials")

        flow_rate(p300, aspirate=20, dispense=20, blow_out=100)
        src = find_materials_well("CPcell", "DW")

        unique_sample = []
        for key in PARAMETERS["Plates"].keys():
            plate=PARAMETERS["Plates"][key]
            if plate["type"] == "TF":
                unique_sample += list(plate["data"].values())
        unique_sample = list(set(unique_sample))
        if "" in unique_sample:
            unique_sample.remove("")

        dest = [
            find_materials_well(name, "DNA", right_well=True)
            for name in unique_sample
        ]

        CP_cell_volume = 40
        ## Mix CP cell
        p300.pick_up_tip()
        for _ in range(2):
            p300.aspirate(25, src)
            p300.dispense(25, src.bottom(z=10))

        # Transfer CP cell to right well of assembly product
        p300.distribute(CP_cell_volume, src.bottom(z=3), dest,
                        new_tip="never", touch_tip=False, disposal_volume=10,
                        blow_out=True, blowout_location="source well",
                        trash=not debug)
        p300.drop_tip()

        # Transfer Assembly Mix to distributed CP cell
        src = [
            find_materials_well(name, "DNA", right_well=False)
            for name in unique_sample
        ]
        reaction_mix_vol = 5
        p20.transfer(reaction_mix_vol, src, dest,
                      new_tip='always',
                      blow_out=False, trash=not debug)

        tc_mod.close_lid()
        protocol.delay(minutes=10)

        ## Thermocycler
        tc_mod.set_block_temperature(
            temperature=42,
            hold_time_seconds=90,
            block_max_volume=reaction_mix_vol + CP_cell_volume
        )
        tc_mod.set_block_temperature(
            temperature = 8,
            hold_time_minutes=2,
            block_max_volume=final_volume
            )
        tc_mod.open_lid()

        src = find_materials_well("LB", "DW").bottom(z=3)

        # Add media for recovery
        #start_time = time.time()
        p300.transfer(60, src, dest,
                    new_tip="always", touch_tip=False, disposal_volume=5,
                    blow_out=False, trash=not debug)
        #protocol.delay(seconds=30)
        #end_time = time.time()

        # Duration time in 8 degree
        #rest_time = 120 - (end_time - start_time)
        #if rest_time > 0:
        #    if rest_time < 30:
        #        protocol.delay(seconds=30)
        #    else:
        #        protocol.delay(seconds=rest_time)

        # Recovery
        tc_mod.set_block_temperature(temperature=37,
                                     hold_time_minutes=int(int(PARAMETERS["Parameters"]["TF_recovery_time"])/2))

        for dest_well in dest:
            p300.pick_up_tip()
            for _ in range(2):
                p300.aspirate(40, dest_well)
                p300.dispense(40, dest_well.bottom(z=5))
            p300.drop_tip()

        tc_mod.set_block_temperature(temperature=37,
                                        hold_time_minutes=int(int(PARAMETERS["Parameters"]["TF_recovery_time"])/2))

        #total_duration = int(PARAMETERS["Parameters"]["TF_recovery_time"]) * 60
        #start_time = time.time()
        #end_time = start_time + total_duration
        #interval = 5 * 60

        #flow_rate(p300, aspirate=40, dispense=40, blow_out=40)
        #if not protocol.is_simulating(): # Avoid simulation error
        #    while time.time() < end_time - 120:
        #        for dest_well in dest:
        #            if time.time() < end_time - 60:
        #                break
        #            p300.pick_up_tip()
        #            for _ in range(2):
        #                p300.aspirate(40, dest_well)
        #                p300.dispense(40, dest_well.bottom(z=5))
        #            p300.drop_tip()

        #        protocol.delay(seconds=interval - len(dest)*10)

        # Spotting
        spotting_volume = 4
        flow_rate(p20, aspirate=8, dispense=15, blow_out=15)
        for i in PARAMETERS["Plates"].keys():
            plate = PARAMETERS["Plates"][i]
            if plate["type"] != "TF":
                continue
            else:
                unique_sample = list(set(plate["data"].values()))
                if "" in unique_sample:
                    unique_sample.remove("")
                for sample in unique_sample:
                    src = find_materials_well(sample, "DNA", right_well=True)
                    dest = [plate['Deck'][well] for well, value in plate["data"].items() if value == sample]
                    p20.pick_up_tip()
                    # Mix Sample
                    for _ in range(3):
                        p20.aspirate(20, src)
                        p20.dispense(20, src.bottom(z=4))
                    # Dispense Sample
                    spotting_dispense(p20, src, dest, spotting_volume)
                    # p20.distribute(spotting_volume, src, dest, new_tip="never", touch_tip=False, trash=not debug)
                    p20.drop_tip()
        tc_mod.deactivate()

    discord_message(f"Protocol End: {time.strftime('%Y-%m-%d %H:%M:%S')}")
