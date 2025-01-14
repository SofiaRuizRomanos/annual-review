def group_c_not_agg(cons, D):

    # 1. find the last element's row and column
    # 2. find the rest of the mass of peope included in that element to reach our %
    # 3. find what share of the mass of people in the cell rest represents
    # 4. create share matrix and multiply by c

    if D.shape != cons.shape:
        raise Exception("Both distributions do not have the same size.")

    locate_matrices = []

    for percentage in [0.05, 0.25]:
        for top in [True, False]:
            beta_matrices = []

            for beta_state in [
                1,
                2,
            ]:  # taking 50% of our percentage for bhigh and 50% from blow
                if beta_state == 1:
                    D_beta = D[:33, :]
                else:
                    D_beta = D[33:, :]

                last_row = D_beta.shape[0]
                last_column = D_beta.shape[1]

                cum_curr_cell = 0
                cum_prev_cell = 0
                last_cell = np.array([last_row, last_column])
                rest_share = 0
                condition_met = False

                if top == True:
                    for row in range(
                        last_row, -1, -1
                    ):  # Iterate over rows from bottom to top
                        for col in range(
                            last_column, -1, -1
                        ):  # Iterate over columns from right to left
                            cum_curr_cell = cum_curr_cell + D_beta[row - 1, col - 1]
                            if (
                                cum_curr_cell == percentage
                            ):  # if in this cell we meet the exact mass percentage we are looking for
                                # stop
                                condition_met = True
                                last_cell = np.array([row - 1, col - 1])
                                rest_share = 1
                                break
                            elif (
                                cum_curr_cell > percentage
                            ):  # if at this cell we meet more than the exact mass percentage we are looking for
                                condition_met = True
                                last_cell = np.array([row - 1, col - 1])
                                rest_share = (percentage - cum_prev_cell) / D_beta[
                                    row - 1, col - 1
                                ]
                                print(rest_share)
                            elif (
                                cum_curr_cell < percentage
                            ):  # if at this cell we meet less than the mass percentage we are looking for
                                # sum and continue
                                cum_prev_cell = cum_curr_cell
                        if condition_met:
                            break
                else:
                    for row in range(
                        0, last_row
                    ):  # Iterate over rows from top to bottom
                        for col in range(
                            0, last_column
                        ):  # Iterate over columns from left to right
                            cum_curr_cell = cum_curr_cell + D_beta[row, col]
                            if (
                                cum_curr_cell == percentage
                            ):  # if in this cell we meet the exact mass percentage we are looking for
                                # stop
                                condition_met = True
                                last_cell = np.array([row, col])
                                rest_share = 1
                                break
                            elif (
                                cum_curr_cell > percentage
                            ):  # if at this cell we meet more than the exact mass percentage we are looking for
                                condition_met = True
                                last_cell = np.array([row, col])
                                # print(percentage)
                                # print(cum_prev_cell)
                                # print(D_beta[row,col])
                                rest_share = (percentage - cum_prev_cell) / D_beta[
                                    row, col
                                ]
                                print("red")
                            elif (
                                cum_curr_cell < percentage
                            ):  # if at this cell we meet less than the mass percentage we are looking for
                                # sum and continue
                                cum_prev_cell = cum_curr_cell
                        if condition_met:
                            break

                provisional_matrix = np.zeros_like(D_beta)
                provisional_matrix[tuple(last_cell)] = rest_share
                print(rest_share)
                # print(provisional_matrix[tuple(last_cell)])
                for r in range(0, last_row):
                    for c in range(0, last_column):
                        if (r > last_cell[0]) or (
                            r == last_cell[0] and c > last_cell[1]
                        ):
                            provisional_matrix[r, c] = 1
                beta_matrices.append(provisional_matrix)
            locate_matrices.append(np.vstack(beta_matrices))

    loc_test = locate_matrices[0]
    c_t10 = np.multiply(cons, locate_matrices[0])
    c_b10 = cons * np.atleast_2d(locate_matrices[1])
    c_t50 = cons * np.atleast_2d(locate_matrices[2])
    c_b50 = cons * np.atleast_2d(locate_matrices[3])
    return loc_test, c_t10, c_b10, c_t50, c_b50


def top_bottom_consumption_2_scrap(c, D):
    if D.shape != c.shape:
        raise Exception("Both distributions do not have the same size.")

    fct_output = []

    for percentage in [0.05, 0.25]:
        for top in [True, False]:
            if top == False:
                percentage = 1 - percentage

            for beta_state in [1, 2]:
                beta_sums = []

                if beta_state == 1:
                    D_beta = D[:33, :]
                    c_beta = c[:33, :]
                else:
                    D_beta = D[33:, :]
                    c_beta = c[33:, :]

                last_row = D_beta.shape[0]
                last_column = D_beta.shape[1]

                pre_sum = 0
                element_percentage = 0
                rest = 0
                last_element = np.array([last_row - 1, last_column - 1])
                condition_met = False
                var_aggregated = 0

                for row in range(
                    last_row - 1, -1, -1
                ):  # Iterate over rows from bottom to top
                    for col in range(
                        last_column - 1, -1, -1
                    ):  # Iterate over columns from right to left
                        element_percentage = pre_sum + D_beta[row, col]
                        if element_percentage == percentage:
                            # if at this cell we meet the exact mass percentage we are looking for
                            var_aggregated += D_beta[row, col] * c_beta[row, col]
                            last_element = np.array([row, col])
                            pre_sum = element_percentage
                            condition_met = True
                            beta_sums.append(var_aggregated)
                            break

                        elif element_percentage > percentage:
                            # if at this cell we meet more than the exact mass percentage we are looking for
                            rest = percentage - pre_sum
                            var_aggregated += rest * c_beta[row, col]
                            last_element = np.array([row, col])
                            condition_met = True
                            beta_sums.append(var_aggregated)
                            break

                        elif element_percentage < percentage:
                            # if at this cell we meet less than the mass percentage we are looking for
                            # sum and continue
                            pre_sum = element_percentage
                            var_aggregated += D_beta[row, col] * c_beta[row, col]
                            last_element = np.array([row, col])

                    if condition_met:
                        break

                beta_sum_total = sum(beta_sums)

            fct_output.append(beta_sum_total)

    c_t10 = fct_output[0]
    c_b10 = fct_output[1]
    c_t50 = fct_output[2]
    c_b50 = fct_output[3]

    return c_t10, c_b10, c_t50, c_b50


def top_bottom_c_agg(c, D):
    if D.shape != c.shape:
        raise Exception("Both distributions do not have the same size.")

    fct_output = []

    for percentage in [0.05, 0.25]:
        for top in [True, False]:
            if top == False:
                percentage = 1 - percentage

            beta_sums = []
            for beta_state in [1, 2]:
                if beta_state == 1:  # s'assurer de bhi et blo
                    D_beta = D[:33, :]
                    c_beta = c[:33, :]
                else:
                    D_beta = D[33:, :]
                    c_beta = c[33:, :]

                last_row = D_beta.shape[0]
                last_column = D_beta.shape[1]

                pre_sum = 0
                element_percentage = 0
                rest = 0
                last_element = np.array([last_row - 1, last_column - 1])
                condition_met = False
                var_aggregated = 0

                for row in range(
                    last_row, -1, -1
                ):  # Iterate over rows from bottom to top
                    for col in range(
                        last_column, -1, -1
                    ):  # Iterate over columns from right to left
                        element_percentage = pre_sum + D_beta[row - 1, col - 1]

                        if element_percentage == percentage:
                            # if at this cell we meet the exact mass percentage we are looking for
                            var_aggregated += (
                                D_beta[row - 1, col - 1] * c_beta[row - 1, col - 1]
                            )
                            last_element = np.array([row - 1, col - 1])
                            pre_sum = element_percentage
                            condition_met = True
                            if top == False:
                                var_aggregated = 1 - var_aggregated
                            beta_sums.append(var_aggregated)
                            break

                        elif element_percentage > percentage:
                            # if at this cell we meet more than the exact mass percentage we are looking for
                            rest = percentage - pre_sum
                            var_aggregated += rest * c_beta[row - 1, col - 1]
                            last_element = np.array([row - 1, col - 1])
                            condition_met = True
                            if top == False:
                                var_aggregated = 1 - var_aggregated
                            beta_sums.append(var_aggregated)
                            break

                        elif element_percentage < percentage:
                            # if at this cell we meet less than the mass percentage we are looking for
                            # sum and continue
                            pre_sum = element_percentage
                            var_aggregated += (
                                D_beta[row - 1, col - 1] * c_beta[row - 1, col - 1]
                            )
                            last_element = np.array([row - 1, col - 1])

                    if condition_met:
                        break

            beta_sum_total = sum(beta_sums)
            fct_output.append(beta_sum_total)

    print(fct_output)
    C_T10 = fct_output[0]
    C_B10 = fct_output[1]
    C_T50 = fct_output[2]
    C_B50 = fct_output[3]

    return C_T10, C_B10, C_T50, C_B50


def top_c_agg(c, D):
    if D.shape != c.shape:
        raise Exception("Both distributions do not have the same size.")

    fct_output = []
    for percentage in [0.05, 0.25]:

        beta_sums = []
        for beta_state in [1, 2]:
            if beta_state == 1:  # s'assurer de bhi et blo
                D_beta = D[:33, :]
                c_beta = c[:33, :]
            else:
                D_beta = D[33:, :]
                c_beta = c[33:, :]

            last_row = D_beta.shape[0]
            last_column = D_beta.shape[1]
            pre_sum = 0
            element_percentage = 0
            rest = 0
            last_element = np.array([last_row - 1, last_column - 1])
            condition_met = False
            var_aggregated = 0

            for row in range(last_row, -1, -1):  # Iterate over rows from bottom to top
                for col in range(
                    last_column, -1, -1
                ):  # Iterate over columns from right to left
                    element_percentage = pre_sum + D_beta[row - 1, col - 1]

                    if element_percentage == percentage:
                        # if at this cell we meet the exact mass percentage we are looking for
                        var_aggregated += (
                            D_beta[row - 1, col - 1] * c_beta[row - 1, col - 1]
                        )
                        last_element = np.array([row - 1, col - 1])
                        pre_sum = element_percentage
                        condition_met = True
                        beta_sums.append(var_aggregated)
                        break

                    elif element_percentage > percentage:
                        # if at this cell we meet more than the exact mass percentage we are looking for
                        rest = percentage - pre_sum
                        var_aggregated += rest * c_beta[row - 1, col - 1]
                        last_element = np.array([row - 1, col - 1])
                        condition_met = True
                        beta_sums.append(var_aggregated)
                        break

                    elif element_percentage < percentage:
                        # if at this cell we meet less than the mass percentage we are looking for
                        # sum and continue
                        pre_sum = element_percentage
                        var_aggregated += (
                            D_beta[row - 1, col - 1] * c_beta[row - 1, col - 1]
                        )
                        last_element = np.array([row - 1, col - 1])

                if condition_met:
                    break
        beta_sum_total = sum(beta_sums)
        fct_output.append(beta_sum_total)

    print(fct_output)
    C_T10 = fct_output[0]
    C_T50 = fct_output[1]

    return C_T10, C_T50


def bottom_c_agg(c, D):
    if D.shape != c.shape:
        raise Exception("Both distributions do not have the same size.")

    fct_output = []
    for percentage in [0.05, 0.25]:

        beta_sums = []
        for beta_state in [1, 2]:
            if beta_state == 1:  # s'assurer de bhi et blo
                D_beta = D[:33, :]
                c_beta = c[:33, :]
            else:
                D_beta = D[33:, :]
                c_beta = c[33:, :]

            last_row = D_beta.shape[0]
            last_column = D_beta.shape[1]
            pre_sum = 0
            element_percentage = 0
            rest = 0
            last_element = np.array([last_row - 1, last_column - 1])
            condition_met = False
            var_aggregated = 0

            for row in range(0, last_row):
                for col in range(0, last_column):
                    element_percentage = pre_sum + D_beta[row, col]

                    if element_percentage == percentage:
                        # if at this cell we meet the exact mass percentage we are looking for
                        var_aggregated += D_beta[row, col] * c_beta[row, col]
                        last_element = np.array([row, col])
                        pre_sum = element_percentage
                        condition_met = True
                        beta_sums.append(var_aggregated)
                        break

                    elif element_percentage > percentage:
                        # if at this cell we meet more than the exact mass percentage we are looking for
                        rest = percentage - pre_sum
                        var_aggregated += rest * c_beta[row, col]
                        last_element = np.array([row, col])
                        condition_met = True
                        beta_sums.append(var_aggregated)
                        break

                    elif element_percentage < percentage:
                        # if at this cell we meet less than the mass percentage we are looking for
                        # sum and continue
                        pre_sum = element_percentage
                        var_aggregated += D_beta[row, col] * c_beta[row, col]
                        last_element = np.array([row, col])

                if condition_met:
                    break
        beta_sum_total = sum(beta_sums)
        fct_output.append(beta_sum_total)

    print(fct_output)
    C_B10 = fct_output[0]
    C_B50 = fct_output[1]

    return C_B10, C_B50
